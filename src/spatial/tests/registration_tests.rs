use std::sync::Arc;

use crate::collector::Count;
use crate::directory::RamDirectory;
use crate::indexer::segment_updater::merge_filtered_segments;
use crate::indexer::{SegmentWriter, SingleSegmentIndexWriter};
use crate::plugin::{PluginMergeContext, PluginWriter, PluginWriterContext, SegmentPlugin};
use crate::query::SpatialQuery;
use crate::schema::document::ErasedValue;
use crate::schema::{Field, FieldEntry, Schema, SpatialOptions, PLANE, SPHERE, STORED};
use crate::spatial::geometry_set::GeometrySet;
use crate::spatial::plane::Plane;
use crate::spatial::plugin::SpatialIndexCatalog;
use crate::spatial::sphere::Sphere;
use crate::spatial::surface::Surface;
use crate::spatial::{Cons, NativeSpatial, Nil, S2Plugin, SpatialInput, SurfaceIndex};
use crate::{DocId, Index, IndexSettings, TantivyDocument, TantivyError};

type Both = Cons<SurfaceIndex<Sphere>, Cons<SurfaceIndex<Plane>, Nil>>;
type Duplicate = Cons<SurfaceIndex<Sphere>, Cons<SurfaceIndex<Sphere>, Nil>>;

fn plugin<Catalog: SpatialIndexCatalog<NativeSpatial>>() -> Arc<dyn SegmentPlugin> {
    Arc::new(S2Plugin::new::<Catalog, NativeSpatial>())
}

fn schema(options: SpatialOptions) -> (Schema, Field) {
    let mut builder = Schema::builder();
    let field = builder.add_spatial_field("geometry", options);
    (builder.build(), field)
}

fn assert_error<T>(result: crate::Result<T>, expected: &str) {
    let error = result.err().expect("registration should be rejected");
    assert!(error.to_string().contains(expected), "{error}");
}

struct ExtensionOwner;

impl SegmentPlugin for ExtensionOwner {
    fn extensions(&self) -> &[&str] {
        &["cells", "edges", "docids"]
    }

    fn create_writer(&self, _ctx: &PluginWriterContext) -> crate::Result<Box<dyn PluginWriter>> {
        unreachable!("extension ownership alone cannot satisfy a spatial field")
    }

    fn merge(&self, _ctx: PluginMergeContext) -> crate::Result<()> {
        unreachable!("extension ownership alone cannot satisfy a spatial field")
    }
}

struct DeclinesSpatial;

impl SpatialInput for DeclinesSpatial {
    fn index_name(_entry: &FieldEntry) -> crate::Result<Option<&str>> {
        Ok(None)
    }

    fn convert<S: Surface>(
        _doc_id: DocId,
        _value: ErasedValue<'_>,
    ) -> crate::Result<Option<GeometrySet<S>>> {
        unreachable!("an unclaimed native spatial field must be rejected at admission")
    }
}

fn invalid_registrations() -> Vec<(Vec<Arc<dyn SegmentPlugin>>, &'static str)> {
    vec![
        (vec![], "requires a registered spatial segment plugin"),
        (vec![plugin::<Nil>()], "'sphere' is not registered"),
        (
            vec![plugin::<Duplicate>()],
            "'sphere' is registered more than once",
        ),
        (vec![plugin::<Both>(), plugin::<Both>()], "more than one"),
        (
            vec![Arc::new(ExtensionOwner)],
            "requires a registered spatial segment plugin",
        ),
        (
            vec![Arc::new(S2Plugin::new::<Both, DeclinesSpatial>())],
            "not supported by the registered plugin input",
        ),
    ]
}

#[test]
fn test_registration_before_metadata_creation() -> crate::Result<()> {
    let (schema, _) = schema(SPHERE);
    for (plugins, expected) in invalid_registrations() {
        let expected = if plugins.len() > 1 {
            "Conflicting segment plugins"
        } else {
            expected
        };
        let directory = RamDirectory::default();
        let mut builder = Index::builder().schema(schema.clone());
        for plugin in plugins {
            builder = builder.register_plugin(plugin);
        }
        assert_error(builder.create(directory.clone()), expected);
        assert!(!Index::exists(&directory)?);
    }
    let directory = RamDirectory::default();
    assert_error(
        Index::create(directory.clone(), schema, IndexSettings::default()),
        "requires a registered spatial segment plugin",
    );
    assert!(!Index::exists(&directory)?);
    Ok(())
}

#[test]
fn test_registration_requires_explicit_name() -> crate::Result<()> {
    let missing: SpatialOptions = serde_json::from_str(r#"{"stored":true}"#)?;
    for options in [SpatialOptions::default(), STORED.into(), missing] {
        let (schema, _) = schema(options);
        assert_error(
            Index::builder()
                .schema(schema)
                .register_plugin(plugin::<Both>())
                .create_in_ram(),
            "requires an explicit spatial index name",
        );
    }
    Ok(())
}

#[test]
fn test_registration_on_existing_builder_and_query() -> crate::Result<()> {
    let (schema, field) = schema(SPHERE);
    let directory = RamDirectory::default();
    Index::builder()
        .schema(schema.clone())
        .register_plugin(plugin::<Both>())
        .create(directory.clone())?;
    let query = SpatialQuery::intersects_bounds(field, [[0.0, 0.0], [1.0, 1.0]]);
    for (plugins, expected) in invalid_registrations() {
        let admission_error = if plugins.len() > 1 {
            "Conflicting segment plugins"
        } else {
            expected
        };
        let mut builder = Index::builder().schema(schema.clone());
        let mut reopened = Index::open(directory.clone())?;
        for plugin in plugins {
            builder = builder.register_plugin(plugin.clone());
            reopened.register_plugin(plugin);
        }
        assert_error(builder.open_or_create(directory.clone()), admission_error);
        // File readers remain available without executable spatial registration.
        let reader = reopened.reader()?;
        assert_error(reader.searcher().search(&query, &Count), expected);
    }
    let reopened = Index::builder()
        .schema(schema)
        .register_plugin(plugin::<Both>())
        .open_or_create(directory)?;
    assert_eq!(reopened.reader()?.searcher().search(&query, &Count)?, 0);
    Ok(())
}

#[test]
fn test_registration_at_write_and_merge_admission() -> crate::Result<()> {
    let (schema, _) = schema(SPHERE);
    let directory = RamDirectory::default();
    let index = Index::builder()
        .schema(schema)
        .register_plugin(plugin::<Both>())
        .create(directory.clone())?;
    // Keep an emitted zero-document source so merge cannot skip admission.
    let writer: SingleSegmentIndexWriter = SingleSegmentIndexWriter::new(index, 50_000_000)?;
    writer.finalize()?;
    for (plugins, expected) in invalid_registrations() {
        let field_error = expected;
        let expected = if plugins.is_empty() {
            "Missing required segment plugin"
        } else if plugins.len() > 1 {
            "Conflicting segment plugins"
        } else {
            expected
        };
        let mut reopened = Index::open(directory.clone())?;
        for plugin in plugins {
            reopened.register_plugin(plugin);
        }
        assert_error(
            reopened.writer_with_num_threads::<TantivyDocument>(1, 50_000_000),
            expected,
        );
        assert_error(
            SingleSegmentIndexWriter::<TantivyDocument>::new(reopened.clone(), 50_000_000),
            expected,
        );
        assert_error(
            SegmentWriter::for_segment(50_000_000, reopened.new_segment()),
            field_error,
        );
        let segments = reopened.searchable_segments()?;
        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].meta().num_docs(), 0);
        let target = RamDirectory::default();
        assert_error(
            merge_filtered_segments(
                &segments,
                IndexSettings::default(),
                vec![None],
                target.clone(),
            ),
            expected,
        );
        assert!(!Index::exists(&target)?);
    }
    Ok(())
}

#[test]
fn test_registration_is_captured_by_readers_and_segments() -> crate::Result<()> {
    let (schema, field) = schema(SPHERE | STORED);
    let directory = RamDirectory::default();
    let mut writer: SingleSegmentIndexWriter = Index::builder()
        .schema(schema.clone())
        .register_plugin(plugin::<Both>())
        .single_segment_index_writer(directory.clone(), 50_000_000)?;
    let document = TantivyDocument::parse_json(
        &schema,
        r#"{"geometry":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]}}"#,
    )?;
    writer.add_document(document)?;
    writer.finalize()?;
    let mut reopened = Index::open(directory)?;
    let before_attachment = reopened.reader_builder();
    let segment = reopened.new_segment();
    reopened.register_plugin(plugin::<Both>());
    let reader = before_attachment.try_into()?;
    let searcher = reader.searcher();
    let stored: TantivyDocument = searcher.doc(crate::DocAddress::new(0, 0))?;
    assert_eq!(stored.get_all(field).count(), 1);
    let query = SpatialQuery::intersects_bounds(field, [[-1.0, -1.0], [2.0, 2.0]]);
    assert_error(
        searcher.search(&query, &Count),
        "requires a registered spatial segment plugin",
    );
    assert_error(
        SegmentWriter::for_segment(50_000_000, segment),
        "requires a registered spatial segment plugin",
    );
    assert_eq!(reopened.reader()?.searcher().search(&query, &Count)?, 1);
    Ok(())
}

#[test]
fn test_registration_routes_mixed_fields() -> crate::Result<()> {
    let mut builder = Schema::builder();
    let sphere = builder.add_spatial_field("sphere", SPHERE);
    let plane = builder.add_spatial_field("plane", PLANE);
    let schema = builder.build();
    let mut writer: SingleSegmentIndexWriter = Index::builder()
        .schema(schema.clone())
        .register_plugin(plugin::<Both>())
        .single_segment_index_writer(RamDirectory::default(), 50_000_000)?;
    writer.add_document(TantivyDocument::parse_json(&schema,
        r#"{"sphere":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},"plane":{"type":"Polygon","coordinates":[[[10,10],[11,10],[11,11],[10,11],[10,10]]]}}"#)?)?;
    let index = writer.finalize()?;
    let searcher = index.reader()?.searcher();
    for (field, matching, absent) in [(sphere, 0.0, 10.0), (plane, 10.0, 0.0)] {
        for (coordinate, expected) in [(matching, 1), (absent, 0)] {
            let query = SpatialQuery::intersects_bounds(
                field,
                [
                    [coordinate - 1.0, coordinate - 1.0],
                    [coordinate + 2.0, coordinate + 2.0],
                ],
            );
            assert_eq!(searcher.search(&query, &Count)?, expected);
        }
    }
    let distance = SpatialQuery::with_predicate(
        plane,
        vec![[10.0, 10.0]],
        crate::query::SpatialPredicate::DistanceWithin(1.0),
    );
    assert!(
        matches!(searcher.search(&distance, &Count), Err(TantivyError::InvalidArgument(message))
        if message == "positive planar distance is not supported")
    );
    let files = index.load_metas()?.list_segment_files();
    for extension in ["cells", "edges", "docids"] {
        assert_eq!(
            files
                .iter()
                .filter(|path| path.extension().and_then(|ext| ext.to_str()) == Some(extension))
                .count(),
            1
        );
    }
    Ok(())
}

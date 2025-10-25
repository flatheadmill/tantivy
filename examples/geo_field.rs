use tantivy::schema::{Schema, SPATIAL, STORED, TEXT};
use tantivy::{Index, IndexWriter, TantivyDocument, TantivyError};

fn main() -> tantivy::Result<()> {
    let mut schema_builder = Schema::builder();
    schema_builder.add_text_field("name", TEXT | STORED);
    schema_builder.add_spatial_field("geometries", STORED | SPATIAL);
    let schema = schema_builder.build();
    let index = Index::create_in_ram(schema.clone());
    let mut index_writer: IndexWriter = index.writer(50_000_000)?;
    let name = schema.get_field("name").unwrap();
    let geometries = schema.get_field("geometries").unwrap();
    let mut hosmer = TantivyDocument::default();
    let hosmer_geometries = (
        vec![vec![vec![(1, 2), (3, 4)]]], // one multi-polygon
        vec![vec![(1, 2), (3, 4)]],       // one multi-line-sting
        vec![(1, 2), (3, 4)],             // oen multi-point
    );
    Ok(())
}

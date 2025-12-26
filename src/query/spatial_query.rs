//! HUSH

use common::BitSet;

use crate::query::explanation::does_not_match;
use crate::query::{BitSetDocSet, Explanation, Query, Scorer, Weight};
use crate::schema::Field;
use crate::spatial::bvh::{search_contains, search_intersects, search_within, Segment};
use crate::spatial::envelope::{Spatial, SpatialI32, Point};
use crate::spatial::point::GeoPoint;
use crate::{DocId, DocSet, Score, TERMINATED};

#[derive(Clone, Copy, Debug)]
/// HUSH
pub enum SpatialQueryType {
    /// HUSH
    Intersects,
    /// HUSH
    Within,
    /// HUSH
    Contains,
}

#[derive(Clone, Copy, Debug)]
/// HUSH
pub struct SpatialQuery {
    field: Field,
    bounds: [GeoPoint; 2],
    query_type: SpatialQueryType,
}

impl SpatialQuery {
    /// HUSH
    pub fn new(field: Field, bounds: [GeoPoint; 2], query_type: SpatialQueryType) -> Self {
        SpatialQuery {
            field,
            bounds: [bounds[0], bounds[1]],
            query_type,
        }
    }
}

impl Query for SpatialQuery {
    fn weight(
        &self,
        _enable_scoring: super::EnableScoring<'_>,
    ) -> crate::Result<Box<dyn super::Weight>> {
        Ok(Box::new(SpatialWeight::new(
            self.field,
            self.bounds,
            self.query_type,
        )))
    }
}

pub struct SpatialWeight {
    field: Field,
    bounds: [GeoPoint; 2],
    query_type: SpatialQueryType,
}

impl SpatialWeight {
    fn new(field: Field, bounds: [GeoPoint; 2], query_type: SpatialQueryType) -> Self {
        SpatialWeight {
            field,
            bounds,
            query_type,
        }
    }
}

impl Weight for SpatialWeight {
    fn scorer(
        &self,
        reader: &crate::SegmentReader,
        boost: crate::Score,
    ) -> crate::Result<Box<dyn super::Scorer>> {
        let spatial_reader = match reader.spatial_fields().get_field(self.field)? {
            Some(reader) => reader,
            None => {
                let empty_bitset = BitSet::with_max_value(reader.max_doc());
                return Ok(Box::new(SpatialScorer::new(boost, empty_bitset, None)));
            }
        };
        let bvh_tree = Segment::<SpatialI32>::new(spatial_reader.get_bytes());
        match self.query_type {
            SpatialQueryType::Intersects => {
                let mut include = BitSet::with_max_value(reader.max_doc());
                let tr = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[0]);
                let bl = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[1]);
                search_intersects(
                    &bvh_tree,
                    bvh_tree.root_offset,
                    &[ bl.y(), bl.x(), tr.y(), tr.x() ],
                    &mut include,
                )?;
                Ok(Box::new(SpatialScorer::new(boost, include, None)))
            }
            SpatialQueryType::Within => {
                let mut include = BitSet::with_max_value(reader.max_doc());
                // TODO dead code.
                let exclude = BitSet::with_max_value(reader.max_doc());
                let tr = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[0]);
                let bl = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[1]);
                search_within(
                    &bvh_tree,
                    bvh_tree.root_offset,
                    &[ bl.y(), bl.x(), tr.y(), tr.x() ],
                    &mut include,
                )?;
                Ok(Box::new(SpatialScorer::new(boost, include, Some(exclude))))
            }
            SpatialQueryType::Contains => {
                let mut include = BitSet::with_max_value(reader.max_doc());
                // TODO dead code.
                let exclude = BitSet::with_max_value(reader.max_doc());
                let tr = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[0]);
                let bl = <SpatialI32 as Spatial>::Point::from_geo(self.bounds[1]);
                search_contains(
                    &bvh_tree,
                    bvh_tree.root_offset,
                    &[ bl.y(), bl.x(), tr.y(), tr.x() ],
                    &mut include,
                )?;
                Ok(Box::new(SpatialScorer::new(boost, include, Some(exclude))))
            }
        }
    }
    fn explain(
        &self,
        reader: &crate::SegmentReader,
        doc: DocId,
    ) -> crate::Result<super::Explanation> {
        let mut scorer = self.scorer(reader, 1.0)?;
        if scorer.seek(doc) != doc {
            return Err(does_not_match(doc));
        }
        let query_type_desc = match self.query_type {
            SpatialQueryType::Intersects => "SpatialQuery::Intersects",
            SpatialQueryType::Within => "SpatialQuery::Within",
            SpatialQueryType::Contains => "SpatialQuery::Contains",
        };
        let score = scorer.score();
        let mut explanation = Explanation::new(query_type_desc, score);
        explanation.add_context(format!(
            "bounds: [({}, {}), ({}, {})]",
            self.bounds[0].lon, self.bounds[0].lat, self.bounds[1].lon, self.bounds[1].lat,
        ));
        explanation.add_context(format!("field: {:?}", self.field));
        Ok(explanation)
    }
}

struct SpatialScorer {
    include: BitSetDocSet,
    exclude: Option<BitSet>,
    doc_id: DocId,
    score: Score,
}

impl SpatialScorer {
    pub fn new(score: Score, include: BitSet, exclude: Option<BitSet>) -> Self {
        let mut scorer = SpatialScorer {
            include: BitSetDocSet::from(include),
            exclude,
            doc_id: 0,
            score,
        };
        scorer.prime();
        scorer
    }
    fn prime(&mut self) {
        self.doc_id = self.include.doc();
        while self.exclude() {
            self.doc_id = self.include.advance();
        }
    }

    fn exclude(&self) -> bool {
        if self.doc_id == TERMINATED {
            return false;
        }
        match &self.exclude {
            Some(exclude) => exclude.contains(self.doc_id),
            None => false,
        }
    }
}

impl Scorer for SpatialScorer {
    fn score(&mut self) -> Score {
        self.score
    }
}

impl DocSet for SpatialScorer {
    fn advance(&mut self) -> DocId {
        if self.doc_id == TERMINATED {
            return TERMINATED;
        }
        self.doc_id = self.include.advance();
        while self.exclude() {
            self.doc_id = self.include.advance();
        }
        self.doc_id
    }

    fn size_hint(&self) -> u32 {
        match &self.exclude {
            Some(exclude) => self.include.size_hint() - exclude.len() as u32,
            None => self.include.size_hint(),
        }
    }

    fn doc(&self) -> DocId {
        self.doc_id
    }
}

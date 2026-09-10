use std::collections::BTreeMap;
use std::ops::Bound;

use tantivy::collector::{DocSetCollector, TopDocs};
use tantivy::query::{
    AllQuery, BooleanQuery, BoostQuery, EmptyQuery, Occur, Query, QueryParser, QueryParserError,
    RangeQuery, SpatialPredicate, SpatialQuery, TermQuery,
};
use tantivy::schema::{Field, IndexRecordOption, Schema, Value, INDEXED, SPHERE, STORED, STRING};
use tantivy::spatial::executor::{PlanNode, SpatialExecutor, SpatialRelation};
use tantivy::{DocAddress, Index, IndexWriter, Score, Searcher, TantivyDocument, Term};

// One kilometer in radians on a sphere with radius 6371 km.
const ONE_KM: f64 = 1.0 / 6371.0;

struct Fixture {
    searcher: Searcher,
    geometry: Field,
    name: Field,
    kind: Field,
    open: Field,
    area: Field,
}

fn square(lon: f64, lat: f64) -> Vec<[f64; 2]> {
    vec![
        [lon - 0.02, lat - 0.02],
        [lon + 0.02, lat - 0.02],
        [lon + 0.02, lat + 0.02],
        [lon - 0.02, lat + 0.02],
        [lon - 0.02, lat - 0.02],
    ]
}

fn term(field: Field, text: &str) -> Box<dyn Query> {
    Box::new(TermQuery::new(
        Term::from_field_text(field, text),
        IndexRecordOption::Basic,
    ))
}

fn and(queries: Vec<Box<dyn Query>>) -> Box<dyn Query> {
    Box::new(BooleanQuery::new(
        queries
            .into_iter()
            .map(|query| (Occur::Must, query))
            .collect(),
    ))
}

impl Fixture {
    fn new(open_closed_trail: bool) -> tantivy::Result<Self> {
        let mut builder = Schema::builder();
        let geometry = builder.add_spatial_field("geometry", SPHERE);
        let name = builder.add_text_field("name", STORED);
        let kind = builder.add_text_field("kind", STRING);
        let open = builder.add_bool_field("open", INDEXED);
        let area = builder.add_u64_field("area", INDEXED);
        let schema = builder.build();
        let index = Index::create_in_ram(schema.clone());
        let mut writer: IndexWriter = index.writer_with_num_threads(1, 50_000_000)?;
        for (label, kind, area, is_open, lon, lat) in [
            ("p_ok", "park", 100u64, false, 0.00, 0.00),
            ("p_hit", "park", 100, false, 1.00, 0.00),
            ("p_small", "park", 1, false, 1.04, 0.00),
            ("p_far", "park", 100, false, 3.00, 0.00),
            ("t_open", "trail", 0, true, 1.02, 0.01),
            ("t_closed", "trail", 0, open_closed_trail, -0.01, 0.01),
            ("r_open", "road", 0, true, 0.01, 0.02),
            ("n_hit", "road", 100, false, 1.01, 0.02),
            // This park is disjoint from t_open, but within one kilometer.
            ("p_band", "park", 100, false, 1.061, 0.00),
        ] {
            let json = serde_json::json!({
                "name": label,
                "kind": kind,
                "open": is_open,
                "area": area,
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [square(lon, lat)],
                },
            });
            writer.add_document(TantivyDocument::parse_json(&schema, &json.to_string())?)?;
        }
        writer.commit()?;
        let reader = index.reader()?;
        Ok(Self {
            searcher: reader.searcher(),
            geometry,
            name,
            kind,
            open,
            area,
        })
    }

    fn parks(&self) -> Box<dyn Query> {
        term(self.kind, "park")
    }

    fn large(&self) -> Box<dyn Query> {
        Box::new(RangeQuery::new(
            Bound::Included(Term::from_field_u64(self.area, 50)),
            Bound::Unbounded,
        ))
    }

    fn open_features(&self, kind: &str) -> Box<dyn Query> {
        and(vec![
            term(self.kind, kind),
            Box::new(TermQuery::new(
                Term::from_field_bool(self.open, true),
                IndexRecordOption::Basic,
            )),
        ])
    }

    fn join(&self, outer: Box<dyn Query>, relation: SpatialRelation) -> Box<dyn Query> {
        Box::new(SpatialExecutor::new(PlanNode::Join {
            field: self.geometry,
            outer: Box::new(PlanNode::Query(outer)),
            inner: Box::new(PlanNode::Query(self.open_features("trail"))),
            relation,
        }))
    }

    fn literal_intersects(&self) -> SpatialQuery {
        SpatialQuery::with_predicate(
            self.geometry,
            square(1.02, 0.01),
            SpatialPredicate::Intersects,
        )
    }

    fn name(&self, address: DocAddress) -> tantivy::Result<String> {
        let doc: TantivyDocument = self.searcher.doc(address)?;
        Ok(doc
            .get_first(self.name)
            .unwrap()
            .as_str()
            .unwrap()
            .to_owned())
    }

    fn names(&self, query: &dyn Query) -> tantivy::Result<Vec<String>> {
        let mut names = self
            .searcher
            .search(query, &DocSetCollector)?
            .into_iter()
            .map(|address| self.name(address))
            .collect::<tantivy::Result<Vec<_>>>()?;
        names.sort();
        Ok(names)
    }

    fn scores(&self, query: &dyn Query) -> tantivy::Result<BTreeMap<String, Score>> {
        self.searcher
            .search(
                query,
                &TopDocs::with_limit(self.searcher.num_docs() as usize).order_by_score(),
            )?
            .into_iter()
            .map(|(score, address)| Ok((self.name(address)?, score)))
            .collect()
    }

    fn assert_names(&self, query: &dyn Query, expected: &[&str]) -> tantivy::Result<()> {
        let mut expected: Vec<String> = expected.iter().map(|name| (*name).to_owned()).collect();
        expected.sort();
        assert_eq!(self.names(query)?, expected);
        Ok(())
    }
}

#[test]
fn test_join_filters_both_inputs() -> tantivy::Result<()> {
    for open_closed_trail in [false, true] {
        let fixture = Fixture::new(open_closed_trail)?;
        for relation in [SpatialRelation::Near(ONE_KM), SpatialRelation::Intersects] {
            let mut expected = vec!["p_hit"];
            if matches!(relation, SpatialRelation::Near(_)) {
                expected.push("p_band");
            }
            if open_closed_trail {
                expected.push("p_ok");
            }
            let query = fixture.join(and(vec![fixture.parks(), fixture.large()]), relation);
            fixture.assert_names(query.as_ref(), &expected)?;
        }
    }
    Ok(())
}

#[test]
fn test_intersects_join_matches_literal_predicate() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let query = fixture.join(Box::new(AllQuery), SpatialRelation::Intersects);
    let joined = fixture.names(query.as_ref())?;
    assert_eq!(joined, fixture.names(&fixture.literal_intersects())?);
    // Use the literal predicate as the reference for the identical polygon.
    let others: Vec<&str> = joined
        .iter()
        .map(String::as_str)
        .filter(|name| *name != "t_open")
        .collect();
    assert_eq!(others, ["n_hit", "p_hit", "p_small"]);
    Ok(())
}

#[test]
fn test_join_boolean_combinations() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let intersects_self = fixture
        .names(&fixture.literal_intersects())?
        .iter()
        .any(|name| name == "t_open");
    for relation in [SpatialRelation::Near(ONE_KM), SpatialRelation::Intersects] {
        let near = matches!(relation, SpatialRelation::Near(_));
        let mut joined = vec!["n_hit", "p_hit", "p_small"];
        if near {
            joined.push("p_band");
        }
        if near || intersects_self {
            joined.push("t_open");
        }
        let join = fixture.join(Box::new(AllQuery), relation);
        fixture.assert_names(join.as_ref(), &joined)?;

        let union = BooleanQuery::new(vec![
            (Occur::Should, fixture.parks()),
            (Occur::Should, join.box_clone()),
        ]);
        let mut expected = vec!["n_hit", "p_band", "p_far", "p_hit", "p_ok", "p_small"];
        if near || intersects_self {
            expected.push("t_open");
        }
        fixture.assert_names(&union, &expected)?;

        let outside = BooleanQuery::new(vec![
            (Occur::Must, fixture.parks()),
            (Occur::MustNot, join.box_clone()),
        ]);
        let mut expected = vec!["p_far", "p_ok"];
        if !near {
            expected.push("p_band");
        }
        fixture.assert_names(&outside, &expected)?;

        let without_roads = BooleanQuery::new(vec![
            (Occur::MustNot, term(fixture.kind, "road")),
            (Occur::Must, join),
        ]);
        joined.retain(|name| *name != "n_hit");
        fixture.assert_names(&without_roads, &joined)?;

        let large_union = and(vec![fixture.large(), Box::new(union)]);
        fixture.assert_names(
            large_union.as_ref(),
            &["n_hit", "p_band", "p_far", "p_hit", "p_ok"],
        )?;
    }
    Ok(())
}

#[test]
fn test_nested_join_uses_inner_result() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    for relation in [SpatialRelation::Near(ONE_KM), SpatialRelation::Intersects] {
        let inner = SpatialExecutor::new(PlanNode::Join {
            field: fixture.geometry,
            outer: Box::new(PlanNode::Query(term(fixture.kind, "trail"))),
            inner: Box::new(PlanNode::Query(fixture.open_features("road"))),
            relation: relation.clone(),
        });
        fixture.assert_names(&inner, &["t_closed"])?;
        let query = SpatialExecutor::new(PlanNode::Join {
            field: fixture.geometry,
            outer: Box::new(PlanNode::Query(and(vec![fixture.parks(), fixture.large()]))),
            inner: Box::new(PlanNode::Query(Box::new(inner))),
            relation,
        });
        fixture.assert_names(&query, &["p_ok"])?;
    }
    Ok(())
}

#[test]
fn test_join_with_empty_inner() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    for relation in [SpatialRelation::Near(ONE_KM), SpatialRelation::Intersects] {
        let query = SpatialExecutor::new(PlanNode::Join {
            field: fixture.geometry,
            outer: Box::new(PlanNode::Query(Box::new(AllQuery))),
            inner: Box::new(PlanNode::Query(Box::new(EmptyQuery))),
            relation,
        });
        fixture.assert_names(&query, &[])?;
    }
    Ok(())
}

#[test]
fn test_join_keeps_boolean_score_contributions() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parks = BoostQuery::new(fixture.parks(), 3.0);
    let park_scores = fixture.scores(&parks)?;
    for relation in [SpatialRelation::Near(ONE_KM), SpatialRelation::Intersects] {
        let mut expected = vec!["p_hit", "p_small"];
        if matches!(relation, SpatialRelation::Near(_)) {
            expected.push("p_band");
        }
        let query = and(vec![
            Box::new(parks.clone()),
            fixture.join(Box::new(AllQuery), relation),
        ]);
        fixture.assert_names(query.as_ref(), &expected)?;
        let scores = fixture.scores(query.as_ref())?;
        expected.sort();
        assert_eq!(
            scores.keys().map(String::as_str).collect::<Vec<_>>(),
            expected
        );
        for (name, score) in scores {
            assert!(park_scores[&name] > 0.0);
            assert!((score - (park_scores[&name] + 1.0)).abs() < 1e-6);
        }
    }
    Ok(())
}

#[test]
fn test_parser_builds_spatial_inner_queries() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    for predicate in ["$within(1km, ", "$intersects("] {
        let nested = format!(
            "kind:park AND area:>=50 AND geometry:{predicate}$query(kind:trail AND \
             geometry:{predicate}$query(kind:road AND open:true))))"
        );
        let empty = format!("geometry:{predicate}$query())");
        for (text, expected) in [(nested, vec!["p_ok"]), (empty, Vec::new())] {
            let ast = query_grammar::parse_query(&text).unwrap();
            let from_text = parser.parse_query(&text).unwrap();
            fixture.assert_names(from_text.as_ref(), &expected)?;
            let from_ast = parser.build_query_from_user_input_ast(ast.clone()).unwrap();
            fixture.assert_names(from_ast.as_ref(), &expected)?;
            let (from_text, errors) = parser.parse_query_lenient(&text);
            assert!(errors.is_empty(), "{errors:?}");
            fixture.assert_names(from_text.as_ref(), &expected)?;
            let (from_ast, errors) = parser.build_query_from_user_input_ast_lenient(ast);
            assert!(errors.is_empty(), "{errors:?}");
            fixture.assert_names(from_ast.as_ref(), &expected)?;
        }
    }
    Ok(())
}

fn assert_parser_error(
    fixture: &Fixture,
    parser: &QueryParser,
    text: &str,
    error: QueryParserError,
    expected_names: &[&str],
) -> tantivy::Result<()> {
    let ast = query_grammar::parse_query(text).unwrap();
    assert_eq!(parser.parse_query(text).unwrap_err(), error);
    assert_eq!(
        parser
            .build_query_from_user_input_ast(ast.clone())
            .unwrap_err(),
        error
    );
    let (query, errors) = parser.parse_query_lenient(text);
    assert_eq!(errors.as_slice(), std::slice::from_ref(&error));
    fixture.assert_names(query.as_ref(), expected_names)?;
    let (query, errors) = parser.build_query_from_user_input_ast_lenient(ast);
    assert_eq!(errors.as_slice(), std::slice::from_ref(&error));
    fixture.assert_names(query.as_ref(), expected_names)
}

#[test]
fn test_parser_reports_inner_unknown_field() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    for predicate in ["$within(1km, ", "$intersects("] {
        for inner in [
            "unknown:trail".to_string(),
            format!("kind:trail AND geometry:{predicate}$query(unknown:road))"),
        ] {
            assert_parser_error(
                &fixture,
                &parser,
                &format!("kind:park AND geometry:{predicate}$query({inner}))"),
                QueryParserError::FieldDoesNotExist("unknown".to_string()),
                &[],
            )?;
        }
    }
    Ok(())
}

#[test]
fn test_parser_keeps_recovered_inner_query() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    for (predicate, expected) in [
        ("$within(1km, ", vec!["p_band", "p_hit", "p_ok", "p_small"]),
        ("$intersects(", vec!["p_hit", "p_ok", "p_small"]),
    ] {
        assert_parser_error(
            &fixture,
            &parser,
            &format!("kind:park AND geometry:{predicate}$query(kind:trail AND unknown:x))"),
            QueryParserError::FieldDoesNotExist("unknown".to_string()),
            &expected,
        )?;
    }
    Ok(())
}

#[test]
fn test_parser_validates_negative_inner_root() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    for (predicate, relation) in [
        ("$within(1km, ", SpatialRelation::Near(ONE_KM)),
        ("$intersects(", SpatialRelation::Intersects),
    ] {
        let reference = SpatialExecutor::new(PlanNode::Join {
            field: fixture.geometry,
            outer: Box::new(PlanNode::Query(fixture.parks())),
            inner: Box::new(PlanNode::Query(Box::new(BooleanQuery::new(vec![
                (Occur::Should, Box::new(AllQuery)),
                (Occur::MustNot, term(fixture.kind, "road")),
            ])))),
            relation: relation.clone(),
        });
        let expected = fixture.names(&reference)?;
        assert!(expected.iter().any(|name| name == "p_ok"));
        if matches!(relation, SpatialRelation::Near(_)) {
            // Every park, including p_far, is in the inner input and is near itself.
            fixture.assert_names(&reference, &["p_band", "p_far", "p_hit", "p_ok", "p_small"])?;
        }
        assert_parser_error(
            &fixture,
            &parser,
            &format!("kind:park AND geometry:{predicate}$query(-kind:road))"),
            QueryParserError::AllButQueryForbidden,
            &expected.iter().map(String::as_str).collect::<Vec<_>>(),
        )?;
    }
    Ok(())
}

#[test]
fn test_parser_separates_inner_syntax_and_resolution_errors() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    let text = "kind:park AND geometry:$within(1km, $query(unknown:x AND kind:))";
    let (ast, grammar_errors) = query_grammar::parse_query_lenient(text);
    assert_eq!(
        grammar_errors,
        vec![query_grammar::LenientError {
            pos: 62,
            message: "expected word".to_string(),
        }]
    );
    assert_eq!(
        parser.parse_query(text).unwrap_err(),
        QueryParserError::SyntaxError(text.to_string())
    );
    assert_eq!(
        parser
            .build_query_from_user_input_ast(ast.clone())
            .unwrap_err(),
        QueryParserError::FieldDoesNotExist("unknown".to_string())
    );
    let (query, errors) = parser.parse_query_lenient(text);
    assert_eq!(
        errors,
        vec![
            QueryParserError::SyntaxError("expected word at position 62".to_string()),
            QueryParserError::FieldDoesNotExist("unknown".to_string()),
        ]
    );
    fixture.assert_names(query.as_ref(), &[])?;
    let (query, errors) = parser.build_query_from_user_input_ast_lenient(ast);
    assert_eq!(
        errors,
        vec![QueryParserError::FieldDoesNotExist("unknown".to_string())]
    );
    fixture.assert_names(query.as_ref(), &[])
}

#[test]
fn test_parser_rejects_unsupported_spatial_forms() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    let parks = ["p_band", "p_far", "p_hit", "p_ok", "p_small"];
    for (spatial, message) in [
        (
            "$contains($query(unknown:trail))",
            "$contains joins are not supported.",
        ),
        (
            "$between(1km, 2km, $query(unknown:trail))",
            "$between is not supported.",
        ),
        ("$between(1km, 2km, 1 2)", "$between is not supported."),
        ("$knn(3, 1 2)", "$knn is not supported."),
    ] {
        assert_parser_error(
            &fixture,
            &parser,
            &format!("kind:park AND geometry:{spatial}"),
            QueryParserError::UnsupportedQuery(message.to_string()),
            &parks,
        )?;
    }
    assert_parser_error(
        &fixture,
        &parser,
        "kind:park AND missing:$contains($query(unknown:x))",
        QueryParserError::FieldDoesNotExist("missing".to_string()),
        &parks,
    )?;
    Ok(())
}

#[test]
fn test_parser_validates_programmatic_spatial_forms() -> tantivy::Result<()> {
    use query_grammar::{SpatialPredicateKind, UserInputAst, UserInputLeaf};

    let fixture = Fixture::new(false)?;
    let parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
    for (predicate, coordinates, message) in [
        (
            SpatialPredicateKind::Knn(3),
            Vec::new(),
            "$knn is not supported.",
        ),
        (
            SpatialPredicateKind::Within(ONE_KM.into()),
            vec![(1.0.into(), 2.0.into())],
            "Spatial joins cannot also specify coordinates.",
        ),
    ] {
        let ast = UserInputAst::and(vec![
            query_grammar::parse_query("kind:park").unwrap(),
            UserInputLeaf::Spatial {
                field: Some("geometry".to_string()),
                predicate,
                coordinates,
                inner_query: Some(Box::new(query_grammar::parse_query("kind:trail").unwrap())),
            }
            .into(),
        ]);
        assert_eq!(
            parser
                .build_query_from_user_input_ast(ast.clone())
                .unwrap_err(),
            QueryParserError::UnsupportedQuery(message.to_string())
        );
        let (query, errors) = parser.build_query_from_user_input_ast_lenient(ast);
        assert_eq!(
            errors,
            vec![QueryParserError::UnsupportedQuery(message.to_string())]
        );
        fixture.assert_names(
            query.as_ref(),
            &["p_band", "p_far", "p_hit", "p_ok", "p_small"],
        )?;
    }
    Ok(())
}

#[test]
fn test_parser_applies_fuzzy_settings_to_inner_query() -> tantivy::Result<()> {
    let fixture = Fixture::new(false)?;
    for (predicate, relation, expected) in [
        (
            "$within(1km, ",
            SpatialRelation::Near(ONE_KM),
            vec!["p_band", "p_hit", "p_small"],
        ),
        (
            "$intersects(",
            SpatialRelation::Intersects,
            vec!["p_hit", "p_small"],
        ),
    ] {
        let mut parser = QueryParser::for_index(fixture.searcher.index(), Vec::new());
        let text = format!("kind:park AND geometry:{predicate}$query(kind:trai AND open:true))");
        let ast = query_grammar::parse_query(&text).unwrap();
        let exact = parser.parse_query(&text).unwrap();
        fixture.assert_names(exact.as_ref(), &[])?;
        parser.set_field_fuzzy(fixture.kind, false, 1, true);
        let reference = fixture.join(fixture.parks(), relation);
        fixture.assert_names(reference.as_ref(), &expected)?;
        let parsed = parser.parse_query(&text).unwrap();
        fixture.assert_names(parsed.as_ref(), &expected)?;
        let built = parser.build_query_from_user_input_ast(ast.clone()).unwrap();
        fixture.assert_names(built.as_ref(), &expected)?;
        let (parsed, errors) = parser.parse_query_lenient(&text);
        assert!(errors.is_empty(), "{errors:?}");
        fixture.assert_names(parsed.as_ref(), &expected)?;
        let (built, errors) = parser.build_query_from_user_input_ast_lenient(ast);
        assert!(errors.is_empty(), "{errors:?}");
        fixture.assert_names(built.as_ref(), &expected)?;
    }
    Ok(())
}

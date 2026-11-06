use tantivy::indexer::doc_id_mapping::{DocIdMapping, SegmentDocIdMapping};
use tantivy::{DocAddress, DocId};

fn writer_mapping(mapping: &DocIdMapping) -> &[DocId] {
    mapping.old_to_new_ids()
}

fn merge_mapping(mapping: &SegmentDocIdMapping) -> (Vec<DocAddress>, usize) {
    (
        mapping.iter_old_doc_addrs().collect(),
        mapping.alive_bitsets().len(),
    )
}

#[test]
fn plugin_document_mappings_are_public() {
    let mapping = DocIdMapping::new_permutation(vec![1, 0]).unwrap();
    assert_eq!(writer_mapping(&mapping), &[1, 0]);

    // The merge mapping is constructed by Tantivy and supplied through PluginMergeContext.
    // Referencing its accessors here verifies that an external plugin can consume them.
    let _: fn(&SegmentDocIdMapping) -> (Vec<DocAddress>, usize) = merge_mapping;
}

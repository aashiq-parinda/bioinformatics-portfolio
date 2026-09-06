"""Unit tests for BLAST parsing, hit filtering, and DataFrame conversions."""

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.blast import (
    filter_blast_hits,
    hits_to_dataframe,
    parse_blast_xml_stream,
)

# Standard valid BLAST XML response with parameters and 2 HSP hits
MOCK_BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.14.0+</BlastOutput_version>
  <BlastOutput_reference>Test</BlastOutput_reference>
  <BlastOutput_db>swissprot</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>Test Protein</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>Test Protein</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|12345</Hit_id>
          <Hit_def>sp|P10275|ANDR_HUMAN Androgen receptor [Homo sapiens]</Hit_def>
          <Hit_accession>P10275</Hit_accession>
          <Hit_len>919</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>185.2</Hsp_bit-score>
              <Hsp_score>470</Hsp_score>
              <Hsp_evalue>1.2e-45</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>600</Hsp_hit-from>
              <Hsp_hit-to>699</Hsp_hit-to>
              <Hsp_identity>85</Hsp_identity>
              <Hsp_positive>95</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>MEVQLGLGRVYPRPPSKTYRGAFQNLFQSVREVIQNPGPRHPEAASAAPPGASLLLLQQQ</Hsp_qseq>
              <Hsp_hseq>MEVQLGLGRVYPRPPSKTYRGAFQNLFQSVREVIQNPGPRHPEAASAAPPGASLLLLQQQ</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|67890</Hit_id>
          <Hit_def>sp|Q00001|DISTANT_ORTHOLOG Distant Receptor [Drosophila melanogaster]</Hit_def>
          <Hit_accession>Q00001</Hit_accession>
          <Hit_len>500</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>35.0</Hsp_bit-score>
              <Hsp_score>80</Hsp_score>
              <Hsp_evalue>0.08</Hsp_evalue>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>40</Hsp_query-to>
              <Hsp_hit-from>20</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_identity>7</Hsp_identity>
              <Hsp_positive>12</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>30</Hsp_align-len>
              <Hsp_qseq>TYRGAFQNLFQSVREVIQNPGPRHPEAASA</Hsp_qseq>
              <Hsp_hseq>TYRGAFQNLFQSVREVIQNPGPRHPEAASA</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>2</Statistics_db-num>
          <Statistics_db-len>1419</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def test_parse_blast_xml_stream() -> None:
    stream = io.StringIO(MOCK_BLAST_XML)
    hits = parse_blast_xml_stream(stream, query_length=100)

    assert len(hits) == 2
    top_hit = hits[0]
    assert top_hit.accession == "P10275"
    assert top_hit.organism == "Homo sapiens"
    assert top_hit.identity_pct == 85.0
    assert top_hit.evalue == 1.2e-45
    assert top_hit.query_coverage_pct == 100.0


def test_filter_blast_hits() -> None:
    stream = io.StringIO(MOCK_BLAST_XML)
    raw_hits = parse_blast_xml_stream(stream, query_length=100)

    filtered = filter_blast_hits(raw_hits, min_identity=50.0, max_evalue=1e-5)
    assert len(filtered) == 1
    assert filtered[0].accession == "P10275"


def test_hits_to_dataframe() -> None:
    stream = io.StringIO(MOCK_BLAST_XML)
    hits = parse_blast_xml_stream(stream, query_length=100)
    df = hits_to_dataframe(hits)

    assert df.shape[0] == 2
    assert "accession" in df.columns
    assert "identity_pct" in df.columns
    assert df.iloc[0]["accession"] == "P10275"

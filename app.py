import re
from collections import defaultdict
from enum import Enum, auto
from io import StringIO
from pathlib import Path
from typing import Generator

import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

st.title('AbChain Merger')

uploaded_file = st.file_uploader("Input file", type=["fa", "fasta"], label_visibility="collapsed")
if uploaded_file is None:
    st.markdown('Or paste sequences directly:')
    input_seq = st.text_area("Input sequences (in FASTA format)", height=500, value='')


def validate_regex(regex, exp_groups=1):
    try:
        if not regex.startswith('^'):
            regex = '^' + regex
        if not regex.endswith('$'):
            regex = regex + '$'

        regex = re.compile(regex)
        if exp_groups is not None and exp_groups != regex.groups:
            st.error('Expected one group in regex')
            return None
    except re.error:
        st.error('Not valid regular expression')
        return None
    return regex


col1, col2, col3 = st.columns([1, 1, 1])
with col1:
    heavy_regex_str = col1.text_input("Heavy regex", "(.*)_HC")
    heavy_regex = validate_regex(heavy_regex_str, 1)

with col2:
    light_regex_str = col2.text_input("Light regex", "(.*)_LC")
    light_regex = validate_regex(light_regex_str, 1)

with col3:
    nonpaired = st.radio(
        "Not paired action 👉",
        options=["remove", "keep"],
    )

if not light_regex or not heavy_regex:
    st.stop()

if uploaded_file:
    source = StringIO(uploaded_file.getvalue().decode("utf-8"))
else:
    source = StringIO(input_seq)


class SeqType(Enum):
    LIGHT = auto()
    HEAVY = auto()
    NONE = auto()


def read_fasta_seq(file) -> dict[str, dict[SeqType, SeqRecord]]:
    stems: dict[str, dict[SeqType, SeqRecord]] = defaultdict(dict)
    for record in SeqIO.parse(file, "fasta"):
        heavy_match = heavy_regex.match(record.id)
        light_match = light_regex.match(record.id)
        if light_match and heavy_match:
            raise ValueError(f'Sequence {record.id} matches both conditions')
        if heavy_match:
            if SeqType.HEAVY in stems[heavy_match.group(1)]:
                raise ValueError(f'Heavy sequence already matched for {record.id}')
            stems[heavy_match.group(1)][SeqType.HEAVY] = record
        elif light_match:
            if SeqType.LIGHT in stems[light_match.group(1)]:
                raise ValueError(f'Light sequence already matched for {record.id}')
            stems[light_match.group(1)][SeqType.LIGHT] = record
        else:
            if SeqType.NONE in stems[record.id]:
                raise ValueError(f'Nonmatched sequence already exists for {record.id}')
            stems[record.id][SeqType.NONE] = record

    return stems


def merge_sequences(sequences: dict[str, dict[SeqType, SeqRecord]],
                    return_nonpaired) -> Generator[SeqRecord, None, None]:
    for grp_id, grp_seq in sequences.items():
        # case 1 - L+H only, can save with group_id
        if SeqType.HEAVY in grp_seq and SeqType.LIGHT in grp_seq and SeqType.NONE not in grp_seq:
            yield SeqRecord(seq=grp_seq[SeqType.LIGHT].seq + grp_seq[SeqType.HEAVY].seq, id=grp_id, description='')
            continue
        # case 1 - L+H + non matching
        # cannot save grouped with same name because nonmacthing will be
        if len(grp_seq) == len(SeqType):
            if return_nonpaired:
                yield grp_seq[SeqType.NONE]
            new_name = grp_id + "MERGED"
            yield SeqRecord(seq=grp_seq[SeqType.LIGHT].seq + grp_seq[SeqType.HEAVY].seq, id=new_name, description='')
            continue
        if return_nonpaired:
            for typ in SeqType:
                if typ in grp_seq:
                    yield grp_seq[typ]


raw_records = read_fasta_seq(source)
merged_records = list(merge_sequences(raw_records, return_nonpaired=nonpaired == 'keep'))
rendered_fasta = ''.join(x.format('fasta') for x in merged_records)

if merged_records and len(merged_records) < 100:
    st.text_area(label='Output fasta', value=rendered_fasta, disabled=True, height=500)
if merged_records:
    if uploaded_file:
        download_path = Path(uploaded_file.name)
        download_filename = str(download_path.with_stem(download_path.stem + '_merged'))
    else:
        download_filename = 'merged.fasta'

    st.download_button(label="Download output as fasta",
                       data=rendered_fasta,
                       file_name=download_filename,
                       mime='application/fasta')

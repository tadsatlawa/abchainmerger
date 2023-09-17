import io
from collections import defaultdict
from enum import Enum, auto
from io import StringIO
from pathlib import Path
from typing import Generator, Optional

import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from streamlit_tags import st_tags

st.set_page_config(layout='wide')
st.title('AbChain Merger')

st.markdown('''
Tool for merging separate light and heavy chain to one entry in fasta file.

Matching light+heavy is done entirely on ids from fasta file. Light and heavy chains are matched if both share the same
*base id* (without light/heavy suffix specified below).

Merged sequence id is:
* *base id* only, if no sequence exist with such id in input; so `A_HC` and `A_LC` will be merged to `A`
* *base id* + "_merged", if sequence with *base id* exists input; so `A`, `A_HC` and `A_LC` will be saved
as `A` and `A_merged`)
* *base id* + "_merged1", if sequences with *base id* and *base id* "_merged" exist in input; so `A`, `A_HC`, `A_LC`
and `A_merged` will be saved as `A` and `A_merged` (original sequences) and `A_merged1`)
* if merged1 also exist then saved as merged2 etc

Remarks:
* **Please refer to your organization policy before pasting internal sequences as there is no responsibility assumed**
* Merging is done only by sequence identifiers, no validation if chains are real antibody sequences is performed
(so you can merge any sequences in fact, not only antibodies)
* If there are two sequences with the same id in input file there is an error
* If multiple sequences have the same base if and different prefixes there is an error (there will be ambiguity
of merging); for example `A_HC` and `A_VH` are incorrect.
* As a consequence of above it is impossible to have `A_HC`, `A_VH`, `A_LC`, `A_VK` and merge `A_LC`+`A_HC`
and `A_VH`+`A_VK`.
''')

col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
with col1:
    ending = st.radio('Ending', ['prefix', 'suffix'], horizontal=True, index=1)
with col2:
    light_cases = st_tags(label='Enter options for light chain:', value=['_LC', '_VK', '_VL'], text='Enter to add more')
    if not light_cases:
        st.error('Please provide values')
with col3:
    heavy_cases = st_tags(label='Enter options for heavy chain:', value=['_HC', '_VH'], text='Enter to add more')
    if not heavy_cases:
        st.error('Please provide values')
with col4:
    nonpaired = st.radio("Not paired action", options=["remove", "keep"], index=1, horizontal=True)
if set(light_cases) & set(heavy_cases):
    st.error('Cannot have the same options for light and heavy!')
    st.stop()
if not light_cases or not heavy_cases:
    st.stop()

uploaded_file = st.file_uploader("Input file", type=["fa", "fasta"], label_visibility="collapsed")
if uploaded_file:
    source = StringIO(uploaded_file.getvalue().decode("utf-8"))
else:
    st.markdown('Or paste sequences directly:')
    input_seq = st.text_area("Input sequences (in FASTA format)", height=500, value='')
    source = StringIO(input_seq)


class SeqType(Enum):
    LIGHT = auto()
    HEAVY = auto()
    NONE = auto()


def get_base(seq: str, suffix_prefix: str, values: list[str]) -> Optional[str]:
    if suffix_prefix == 'prefix':
        for val in values:
            if seq.startswith(val):
                return seq[len(val):]
    if suffix_prefix == 'suffix':
        for val in values:
            if seq.endswith(val):
                return seq[:-len(val)]
    return None


def read_fasta_seq(file) -> dict[str, dict[SeqType, SeqRecord]]:
    stems: dict[str, dict[SeqType, SeqRecord]] = defaultdict(dict)
    assert ending
    for record in SeqIO.parse(file, "fasta"):
        heavy_base = get_base(record.id, ending, heavy_cases)
        if heavy_base:
            if SeqType.HEAVY in stems[heavy_base]:
                raise ValueError(f'Heavy sequence already matched for {record.id}')
            stems[heavy_base][SeqType.HEAVY] = record
            continue
        light_base = get_base(record.id, ending, light_cases)
        if light_base:
            if SeqType.LIGHT in stems[light_base]:
                raise ValueError(f'Light sequence already matched for {record.id}')
            stems[light_base][SeqType.LIGHT] = record
            continue
        # neither matched
        if SeqType.NONE in stems[record.id]:
            # if this happens input fasta is wrong and has duplicated ids
            raise ValueError(f'Nonmatched sequence already exists for {record.id}')
        stems[record.id][SeqType.NONE] = record

    return stems


def find_new_name(sequences: dict[str, dict[SeqType, SeqRecord]], current_name: str):
    new_name = current_name + '_merged'
    current_id = 0
    while new_name in sequences:
        current_id += 1
        new_name = current_name + '_merged' + str(current_id)
    return new_name


def merge_sequences(sequences: dict[str, dict[SeqType, SeqRecord]],
                    return_nonpaired) -> Generator[SeqRecord, None, None]:
    for grp_id, grp_seq in sequences.items():
        # case 1 - L+H only, can save with group_id
        if SeqType.HEAVY in grp_seq and SeqType.LIGHT in grp_seq and SeqType.NONE not in grp_seq:
            yield SeqRecord(seq=grp_seq[SeqType.LIGHT].seq + grp_seq[SeqType.HEAVY].seq, id=grp_id, description='')
            continue
        # case 2 - L+H + non matching
        # cannot save grouped with same name because nonmacthing will be
        if len(grp_seq) == len(SeqType):
            if return_nonpaired:
                yield grp_seq[SeqType.NONE]
            new_name = find_new_name(sequences, grp_id)
            yield SeqRecord(seq=grp_seq[SeqType.LIGHT].seq.strip() + grp_seq[SeqType.HEAVY].seq.strip(),
                            id=new_name,
                            description='')
            continue
        if return_nonpaired:
            for typ in SeqType:
                if typ in grp_seq:
                    yield grp_seq[typ]


raw_records = read_fasta_seq(source)
merged_records = list(merge_sequences(raw_records, return_nonpaired=nonpaired == 'keep'))
rendered_fasta = io.StringIO()
SeqIO.write(merged_records, rendered_fasta, 'fasta')

if merged_records:
    if len(merged_records) < 100:
        st.text_area(label='Output fasta',
                     value=rendered_fasta.getvalue(),
                     disabled=True,
                     height=len(merged_records) * 50 + 20)
    else:
        st.write('Too large output to display, please download file.')

    if uploaded_file:
        download_path = Path(uploaded_file.name)
        download_filename = str(download_path.with_stem(download_path.stem + '_merged'))
    else:
        download_filename = 'merged.fasta'

    st.download_button(label="Download output as fasta",
                       data=rendered_fasta.getvalue(),
                       file_name=download_filename,
                       mime='application/fasta')
else:
    st.info('No output generated')

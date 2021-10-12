from collections import defaultdict
from os import getcwd

from vcf import Reader, Writer
from vcf.model import _Record


def get_qual(record1, record2):
    if record1.QUAL is not None:
        if record2.QUAL is not None:
            # both are not None
            if record1.QUAL > record2.QUAL:
                qual = record1.QUAL
            else:
                qual = record2.QUAL
        else:
            # 2nd is None
            qual = record1.QUAL
    else:
        if record2.QUAL is not None:
            # 1st is None
            qual = record2.QUAL
        else:
            # Both are none
            qual = record2.QUAL
    return qual


def get_filter(record1, record2):
    if record1.FILTER is not None:
        # record1 filter is not None
        if record2.FILTER is not None:
            # record2 filter is not None
            if isinstance(record1.FILTER, str) and isinstance(record2.FILTER, str):
                # Both are of type str, so merge and send
                return record1.FILTER + record2.FILTER
            else:
                # Priority to the list one
                if isinstance(record1.FILTER, list):
                    return record1.FILTER
                elif isinstance(record2.FILTER, list):
                    return record2.FILTER
                else:
                    return record1.FILTER
        else:
            return record1.FILTER
    elif record2.FILTER is not None:
        return record2.FILTER
    else:
        return record1.FILTER


def get_samples(record1, record2):
    sample_list = list()
    if len(record1.samples) > 0:
        # sample of rec1 has more than one element
        for sample in record1.samples:
            sample_list.append(sample)
        if len(record2.samples) > 0:
            # sample of rec2 has more than one element
            for sample in record2.samples:
                sample_list.append(sample)
    elif len(record2.samples) > 0:
        # sample of rec2 has more than one element
        for sample in record2.samples:
            sample_list.append(sample)
    else:
        # sample of both records has zero elements
        return record1.samples
    return sample_list


def get_info(record1, record2):
    dd = defaultdict(list)
    common_keys = set.intersection(set(record1.INFO.keys()), set(record2.INFO.keys()))
    for k, v in record1.INFO.items():
        if k in common_keys:
            k = f"Freebayes_{k}"
        dd[k].append(v)
    for k, v in record2.INFO.items():
        if k in common_keys:
            k = f"VarScan_{k}"
        dd[k].append(v)
    dd["calledBy"] = "Freebayes+VarScan"
    return dd


def get_format(record1, record2):
    r1_format = record1.FORMAT.split(":")
    r2_format = record2.FORMAT.split(":")
    common_keys = set.intersection(set(r1_format), set(r2_format))
    format_value_list = list()
    for i in r1_format:
        if i in common_keys:
            format_value_list.append(f"Freebayes_{i}")
        else:
            format_value_list.append(i)
    for i in r2_format:
        if i in common_keys:
            format_value_list.append(f"VarScan_{i}")
        else:
            format_value_list.append(i)
    format_value = ":".join(format_value_list)
    return format_value


def merge_vcf(record1, record2):
    chrom = record1.CHROM  # Since these will contain same chromosomes
    pos = record1.POS  # We're merging for same pos
    identifier = record1.ID
    ref = record1.REF
    alt = record1.ALT
    qual = get_qual(record1, record2)
    filter_value = get_filter(record1, record2)
    info_value = get_info(record1, record2)
    format_value = get_format(record1, record2)

    samples = get_samples(record1, record2)
    record = _Record(
        CHROM=chrom,
        POS=pos,
        ID=identifier,
        REF=ref,
        ALT=alt,
        QUAL=qual,
        FILTER=filter_value,
        INFO=info_value,
        FORMAT=format_value,
        sample_indexes=samples,
    )
    record.samples = samples
    return record


def first_pass(vcf1, vcf2):
    counter1, counter2 = defaultdict(), defaultdict()
    for rec in vcf1:
        if rec.POS not in counter1:
            counter1[rec.POS] = rec
    for rec in vcf2:
        if rec.POS not in counter2:
            counter2[rec.POS] = rec
    return counter1, counter2


def second_pass(counter1, counter2):
    final_output = list()
    common_keys = set.intersection(set(counter1.keys()), set(counter2.keys()))

    for key in common_keys:
        rec_main = merge_vcf(counter1[key], counter2[key])
        final_output.append(rec_main)

    for pos1, rec1 in counter1.items():
        if pos1 not in common_keys:
            rec1.INFO["calledBy"] = "Freebayes"
            final_output.append(rec1)

    for pos2, rec2 in counter2.items():
        if pos2 not in common_keys:
            rec2.INFO["calledBy"] = "VarScan"
            final_output.append(rec2)

    return final_output


def get_header(vcf1, vcf2, output_path):
    vcf1.metadata = {**vcf1.metadata, **vcf2.metadata}
    vcf1.alts = {**vcf1.alts, **vcf2.alts}
    vcf1.contigs = {**vcf1.contigs, **vcf2.contigs}
    vcf1.filename = output_path
    vcf1.filters = {**vcf1.filters, **vcf2.filters}
    vcf1.formats = {**vcf1.formats, **vcf2.formats}
    vcf1.infos = {**vcf1.infos, **vcf2.infos}
    vcf1.samples = vcf1.samples + vcf2.samples


def main():
    fpath1, fpath2 = f"{getcwd()}/freebayes_raw.vcf", f"{getcwd()}/varscan_raw.vcf"
    output_path = f"{getcwd()}/out_final.vcf"
    vcf1 = Reader(open(fpath1))
    vcf2 = Reader(open(fpath2))
    get_header(vcf1, vcf2, output_path)
    counter1, counter2 = first_pass(vcf1, vcf2)
    final_output = second_pass(counter1, counter2)
    vcf_writer = Writer(open(output_path, "w"), vcf1)
    for record in final_output:
        vcf_writer.write_record(record)


main()

import sys
import pysam
import parasail


def hdr_to_shv(hdr_record, open_pen=10, extend_pen=1, match_score=2, mismatch_score=-1):
    ref_hdr, alt_hdr = hdr_record.alleles
    ref_start = hdr_record.pos
    hdr_id = hdr_record.id
    
    aln = parasail.sg_trace_striped_32(
        alt_hdr[1:], ref_hdr[1:],
        open_pen, 
        extend_pen,
        parasail.matrix_create("ACGT", match_score, mismatch_score)
    )
    
    snp_counter = 1
    ins_counter = 1
    del_counter = 1

    ref_pos = ref_start + 1
    ref_idx = 1
    alt_idx = 1
    for iop in aln.cigar.seq:
        ln, op = aln.cigar.decode_len(iop), aln.cigar.decode_op(iop).decode()
        if op == '=':
            # match, no VCF record emitted
            ref_idx += ln
            alt_idx += ln
            ref_pos += ln
        elif op == 'X':
            # mismatch, emit SNP
            for i in range(ln):
                yield ref_pos + i, ref_hdr[ref_idx + i], alt_hdr[alt_idx + i], f'{hdr_id}_SNP{snp_counter}'
                snp_counter += 1
            ref_pos += ln
            ref_idx += ln
            alt_idx += ln
        elif op == 'D':
            yield ref_pos - 1, ref_hdr[ref_idx - 1: ref_idx + ln], alt_hdr[alt_idx - 1], f'{hdr_id}_DEL{del_counter}'
            del_counter += 1
            ref_pos += ln
            ref_idx += ln
        elif op == 'I':
            yield ref_pos - 1, ref_hdr[ref_idx - 1], alt_hdr[alt_idx - 1: alt_idx + ln], f'{hdr_id}_INS{ins_counter}'
            ins_counter += 1
            alt_idx += ln


def convert_syri_output(syri_fn, output_vcf_fn):
    with open(output_vcf_fn, 'w') as o:
        o.write(
            '##fileformat=VCFv4.3\n'
            '##source=syri_for_star\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPARENT2\n'
        )
        with pysam.VariantFile(syri_fn) as v:
            for r in v.fetch():
                if r.info['VarType'] ==  'ShV' and r.ref != 'N' and r.info['Parent'].startswith('SYN'):
                    if not r.id.startswith('HDR'):
                        outrecord = (
                            f'{r.chrom}\t{r.pos}\t{r.id}\t{r.alleles[0]}\t{r.alleles[1]}\t'
                            f'.\tPASS\tMTD=syri\tGT\t0|1\n'
                        )
                        o.write(outrecord)
                    else:
                        for pos, ref, alt, id_ in hdr_to_shv(r):
                            outrecord = (
                                f'{r.chrom}\t{pos}\t{id_}\t{ref}\t{alt}\t'
                                f'.\tPASS\tMTD=syri\tGT\t0|1\n'
                            )
                            o.write(outrecord)


if __name__ == '__main__':
    convert_syri_output(sys.argv[1], sys.argv[2])
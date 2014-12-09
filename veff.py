from __future__ import division, print_function
import collections
import operator


import petl.fluent as etl
import petlx.gff3
import pyfasta
from Bio.Seq import Seq


class Genome(object):

    def __init__(self,
                 fasta_path,
                 gff3_path,
                 seqid=None):
        """
        An annotated reference genome.

        Parameters
        ----------

        fasta_path : string
            Path to reference genome FASTA file.
        gff3_path : string
            Path to genome annotations GFF3 file.

        """

        # store initialisation parameters
        self._fasta_path = fasta_path
        self._gff3_path = gff3_path
        self._seqid = seqid

        # setup access to reference sequence
        self._fasta = pyfasta.Fasta(fasta_path)

        # setup access to GFF3 as a table
        self._tbl_features = (etl
            .fromgff3(gff3_path)
            .unpackdict('attributes', ['ID', 'Parent'])
            .rename({'ID': 'feature_id',
                     'Parent': 'parent_id',
                     'end': 'stop'})
        )

        # limit data to a single chromosome
        if seqid is not None:
            self._tbl_features = self._tbl_features.eq('seqid', seqid)

        # index features by ID
        self._idx_feature_id = self._tbl_features.recordlookupone('feature_id')

        # index features by parent ID
        self._idx_parent_id = self._tbl_features.recordlookup('parent_id')

        # index features by genomic location
        self._idx_location = self._tbl_features.facetintervalrecordlookup(
            'seqid', 'start', 'stop', proximity=1
        )

    def get_feature(self, feature_id):
        return self._idx_feature_id[feature_id]

    def get_children(self, feature_id):
        return self._idx_parent_id[feature_id]

    def find(self, chrom, start, stop):
        return self._idx_location[chrom].find(start, stop)

    def get_ref_allele_coords(self, chrom, pos, ref):

        # N.B., use one-based inclusive coordinate system (like GFF3) throughout
        ref_start = pos
        ref_stop = pos + len(ref) - 1

        # check the reference allele matches the reference sequence
        ref_seq = self._fasta.sequence({'chr': chrom,
                                        'start': ref_start,
                                        'stop': ref_stop})
        ref_seq = str(ref_seq).lower()
        assert ref_seq == ref.lower(), \
            'reference allele does not match reference sequence, ' \
            'expected %r, found %r' % (ref_seq, ref.lower())

        return ref_start, ref_stop


Effect = collections.namedtuple('Effect', ('effect',
                                           'impact',
                                           'chrom',
                                           'pos',
                                           'ref',
                                           'alt',
                                           'ref_start',
                                           'ref_stop',
                                           'gene_id',
                                           'gene_start',
                                           'gene_stop',
                                           'gene_strand',
                                           'transcript_id',
                                           'transcript_start',
                                           'transcript_stop',
                                           'transcript_strand',
                                           'cds_id',
                                           'cds_start',
                                           'cds_stop',
                                           'cds_strand',
                                           'ref_cds_start',
                                           'ref_cds_stop',
                                           'ref_start_phase',
                                           'ref_codon',
                                           'alt_codon',
                                           'codon_change',
                                           'aa_pos',
                                           'ref_aa',
                                           'alt_aa',
                                           'aa_change'))
null_effect = Effect(*([None] * len(Effect._fields)))


def get_effects(genome, chrom, pos, ref, alt,
                gff3_gene_types={'gene', 'pseudogene'},
                gff3_transcript_types={'mRNA', 'pseudogenic_transcript'},
                gff3_cds_types={'CDS', 'pseudogenic_exon'}):
    """TODO

    Parameters
    ----------

    gff3_gene_types : list of strings, optional
        Feature types to consider as genes.
    gff3_transcript_types : list of strings, optional
        Feature types to consider as transcripts.
    gff3_exon_types : list of strings, optional
        Feature types to consider as exons.

    Returns
    -------

    An `Effect` generator.

    """

    # ensure types and case
    ref = str(ref).upper()
    alt = str(alt).upper()

    # obtain start and stop coordinates of the reference allele
    ref_start, ref_stop = genome.get_ref_allele_coords(chrom, pos, ref)

    # setup common effect parameters
    base_effect = null_effect._replace(chrom=chrom,
                                       pos=pos,
                                       ref=ref,
                                       alt=alt,
                                       ref_start=ref_start,
                                       ref_stop=ref_stop)

    # find overlapping genome features
    features = genome.find(chrom, ref_start, ref_stop)

    # filter to find overlapping genes
    genes = [f for f in features if f.type in gff3_gene_types]

    if not genes:
        for effect in _get_intergenic_effects(genome, base_effect):
            yield effect

    else:
        for gene in genes:
            for effect in _get_gene_effects(genome, base_effect, gene,
                                            gff3_transcript_types,
                                            gff3_cds_types):
                yield effect


def _get_intergenic_effects(genome, base_effect):

    # TODO
    # UPSTREAM and DOWNSTREAM

    # the variant is in an intergenic region
    effect = base_effect._replace(effect='INTERGENIC',
                                  impact='MODIFIER')
    yield effect


def _get_gene_effects(genome, base_effect, gene, gff3_transcript_types,
                      gff3_cds_types):

    # setup common effect parameters
    base_effect = base_effect._replace(gene_id=gene.feature_id,
                                       gene_start=gene.start,
                                       gene_stop=gene.stop,
                                       gene_strand=gene.strand)

    # obtain transcripts that are children of the current gene
    transcripts = [t for t in genome.get_children(gene.feature_id)
                   if t.type in gff3_transcript_types]

    if not transcripts:

        # the variant hits a gene, but no transcripts within the gene
        effect = base_effect._replace(effect='INTRAGENIC',
                                      impact='MODIFIER')
        yield effect

    else:

        for transcript in transcripts:
            for effect in _get_transcript_effects(genome, base_effect,
                                                  transcript, gff3_cds_types):
                yield effect


def _get_transcript_effects(genome, base_effect, transcript, gff3_cds_types):

    # setup common effect parameters
    base_effect = base_effect._replace(
        transcript_id=transcript.feature_id,
        transcript_start=transcript.start,
        transcript_stop=transcript.stop,
        transcript_strand=transcript.strand
    )

    # convenience
    ref_start = base_effect.ref_start
    ref_stop = base_effect.ref_stop
    transcript_start = transcript.start
    transcript_stop = transcript.stop

    # compare start and stop of reference allele to start
    # and stop of current transcript

    if ref_stop < transcript_start:

        # TODO
        # variant hits a gene but misses the current transcript, falling
        # upstream
        effect = base_effect._replace(effect='TODO')
        yield effect

    elif ref_start > transcript_stop:

        # TODO
        # variant hits a gene but misses the current transcript, falling
        # downstream
        effect = base_effect._replace(effect='TODO')
        yield effect

    elif ref_start < transcript_start <= ref_stop <= transcript_stop:

        # TODO
        # reference allele overhangs the start of the current transcript
        effect = base_effect._replace(effect='TODO')
        yield effect

    elif transcript_start <= ref_start <= transcript_stop < ref_stop:

        # TODO
        # reference allele overhangs the end of the current transcript
        effect = base_effect._replace(effect='TODO')
        yield effect

    elif ref_start < transcript_start <= transcript_stop < ref_stop:

        # TODO
        # reference allele entirely overlaps the current transcript and
        # overhangs at both ends
        effect = base_effect._replace(effect='TODO')
        yield effect

    else:

        # reference allele falls within current transcript
        assert transcript_start <= ref_start <= ref_stop <= transcript_stop
        yield _get_within_transcript_effect(genome, base_effect, transcript,
                                            gff3_cds_types)


def _get_within_transcript_effect(genome, base_effect, transcript,
                                  gff3_cds_types):

    # convenience
    ref_start = base_effect.ref_start
    ref_stop = base_effect.ref_stop

    # obtain coding sequences that are children of the current transcript
    cdss = [f for f in genome.get_children(transcript.feature_id)
            if f.type in gff3_cds_types]

    if not cdss:

        # TODO
        # the variant hits a transcript, but there are no CDSs within the
        # transcript
        effect = base_effect._replace(effect='TODO')
        return effect

    else:

        # find coding sequence that overlaps the reference allele
        overlapping_cdss = [cds for cds in cdss
                            if cds.start <= ref_start <= cds.stop
                            or cds.start <= ref_stop <= cds.stop
                            or ref_start < cds.start <= cds.stop < ref_stop]

        if not overlapping_cdss:

            # TODO SPLICE_SITE effects

            # variant hits an intron (technically, hits no exon in the
            # transcript)
            effect = base_effect._replace(effect='INTRON',
                                          impact='MODIFIER')
            return effect

        elif len(overlapping_cdss) > 1:

            # TODO
            # variant overlaps more than one exon
            effect = base_effect._replace(effect='TODO')
            return effect

        else:

            # variant overlaps a single exon
            assert len(overlapping_cdss) == 1
            cds = overlapping_cdss[0]

            return _get_cds_effect(genome, base_effect, cds, cdss)


def _get_cds_effect(genome, base_effect, cds, cdss):

    # setup common effect parameters
    base_effect = base_effect._replace(
        cds_id=cds.feature_id,
        cds_start=cds.start,
        cds_stop=cds.stop,
        cds_strand=cds.strand
    )

    # convenience
    ref_start = base_effect.ref_start
    ref_stop = base_effect.ref_stop
    cds_start = cds.start
    cds_stop = cds.stop

    if ref_start < cds_start <= ref_stop <= cds_stop:

        # TODO
        # reference allele overhangs the start of the current exon
        effect = base_effect._replace(effect='TODO')
        return effect

    elif cds_start <= ref_start <= cds_stop < ref_stop:

        # TODO
        # reference allele overhangs the end of the current transcript
        effect = base_effect._replace(effect='TODO')
        return effect

    elif ref_start < cds_start <= cds_stop < ref_stop:

        # TODO
        # reference allele entirely overlaps the current exon and
        # overhangs at both ends
        effect = base_effect._replace(effect='TODO')
        return effect

    else:

        # reference allele falls within current transcript
        assert cds_start <= ref_start <= ref_stop <= cds_stop
        return _get_within_cds_effect(genome, base_effect, cds, cdss)


def _get_within_cds_effect(genome, base_effect, cds, cdss):

    # convenience
    chrom = base_effect.chrom
    pos = base_effect.pos
    ref = base_effect.ref
    alt = base_effect.alt

    # obtain amino acid change
    ref_cds_start, ref_cds_stop, ref_start_phase, ref_codon, alt_codon, \
        aa_pos, ref_aa, alt_aa = _get_aa_change(genome, chrom, pos, ref, alt,
                                                cds, cdss)

    # setup common effect parameters
    base_effect = base_effect._replace(
        ref_cds_start=ref_cds_start,
        ref_cds_stop=ref_cds_stop,
        ref_start_phase=ref_start_phase,
        ref_codon=ref_codon,
        alt_codon=alt_codon,
        codon_change='%s/%s' % (ref_codon, alt_codon),
        aa_pos=aa_pos,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        aa_change='%s%s%s' % (ref_aa, aa_pos, alt_aa)
    )

    if len(ref) == 1 and len(alt) == 1:

        # SNPs

        if ref_aa == alt_aa:

            # TODO SYNONYMOUS_START and SYNONYMOUS_STOP

            # variant causes a codon that produces the same amino acid
            # e.g.: Ttg/Ctg, L/L
            effect = base_effect._replace(effect='SYNONYMOUS_CODING',
                                          impact='LOW')

        elif ref_aa == 'M' and ref_cds_start == 0:

            # variant causes start codon to be mutated into a non-start codon.
            # e.g.: aTg/aGg, M/R
            effect = base_effect._replace(effect='START_LOST',
                                          impact='HIGH')

        elif ref_aa == '*':

            # variant causes stop codon to be mutated into a non-stop codon
            # e.g.: Tga/Cga, */R
            effect = base_effect._replace(effect='STOP_LOST',
                                          impact='HIGH')

        elif alt_aa == '*':

            # variant causes a STOP codon e.g.: Cag/Tag, Q/*
            effect = base_effect._replace(effect='STOP_GAINED',
                                          impact='HIGH')

        else:

            # TODO NON_SYNONYMOUS_START and NON_SYNONYMOUS_STOP

            # variant causes a codon that produces a different amino acid
            # e.g.: Tgg/Cgg, W/R
            effect = base_effect._replace(effect='NON_SYNONYMOUS_CODING',
                                          impact='MODERATE')

    else:

        # INDELs and MNPs

        indel_size = len(ref) - len(alt)

        if indel_size % 3:

            # insertion or deletion causes a frame shift
            # e.g.: An indel size is not multple of 3
            effect = base_effect._replace(effect='FRAME_SHIFT',
                                          impact='HIGH')

        elif indel_size > 0:

            # insertions

            if ref_start_phase > 0:

                # one or many codons are inserted
                # e.g.: An insert multiple of three in a codon boundary
                effect = base_effect._replace(effect='CODON_INSERTION',
                                              impact='MODERATE')

            else:

                # one codon is changed and one or many codons are inserted
                # e.g.: An insert of size multiple of three, not at codon
                # boundary
                effect = base_effect._replace(
                    effect='CODON_CHANGE_PLUS CODON_INSERTION',
                    impact='MODERATE'
                )

        elif indel_size < 0:

            # deletions

            if ref_start_phase > 0:

                # one or many codons are deleted
                # e.g.: A deletions multiple of three in a codon boundary
                effect = base_effect._replace(effect='CODON_DELETION',
                                              impact='MODERATE')

            else:

                # one codon is changed and one or many codons are deleted
                # e.g.: A deletion of size multiple of three, not at codon
                # boundary
                effect = base_effect._replace(
                    effect='CODON_CHANGE_PLUS CODON_DELETION',
                    impact='MODERATE'
                )

        else:

            # MNPs
            assert indel_size == 0
            effect = base_effect._replace(effect='CODON_CHANGE',
                                          impact='MODERATE')

    # TODO check codon changes are non-synonymous?

    return effect


def _get_aa_change(genome, chrom, pos, ref, alt, cds, cdss):

    # obtain codon change
    ref_cds_start, ref_cds_stop, ref_start_phase, ref_codon, alt_codon = \
        _get_codon_change(genome, chrom, pos, ref, alt, cds, cdss)

    # translate codon change to amino acid change
    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())
    aa_pos = (ref_cds_start // 3) + 1

    return ref_cds_start, ref_cds_stop, ref_start_phase, ref_codon, alt_codon, \
        aa_pos, ref_aa, alt_aa


def _get_codon_change(genome, chrom, pos, ref, alt, cds, cdss):

    # convenience
    fasta = genome._fasta

    # obtain reference allele coords relative to coding sequence
    ref_start, ref_stop = genome.get_ref_allele_coords(chrom, pos, ref)
    ref_cds_start, ref_cds_stop = _get_coding_position(ref_start, ref_stop,
                                                       cds, cdss)

    # calculate position of reference allele start within codon
    ref_start_phase = ref_cds_start % 3

    if cds.strand == '+':

        # obtain any previous nucleotides to complete the first codon
        prefix = fasta.sequence({'chr': chrom,
                                 'start': ref_start - ref_start_phase,
                                 'stop': ref_start - 1})
        prefix = str(prefix).lower()

        # begin constructing reference and alternate codon sequences
        ref_codon = prefix + ref
        alt_codon = prefix + alt

        # obtain any subsequence nucleotides to complete the last codon
        if len(ref_codon) % 3:
            ref_stop_phase = len(ref_codon) % 3
            suffix = fasta.sequence({'chr': chrom,
                                     'start': ref_stop + 1,
                                     'stop': ref_stop + 3 - ref_stop_phase})
            suffix = str(suffix).lower()
            ref_codon += suffix

        if len(alt_codon) % 3:
            alt_stop_phase = len(alt_codon) % 3
            suffix = fasta.sequence({'chr': chrom,
                                     'start': ref_stop + 1,
                                     'stop': ref_stop + 3 - alt_stop_phase})
            suffix = str(suffix).lower()
            alt_codon += suffix

    else:

        # N.B., we are on the reverse strand, so position reported for
        # variant is actually position at the *end* of the reference allele
        # which is particularly important for deletions

        # we will construct everything for the forward strand (i.e., back-to-
        # front) then take reverse complement afterwards at the end of this
        # code block

        # obtain any previous nucleotides to complete the first codon
        prefix = fasta.sequence({'chr': chrom,
                                 'start': ref_stop + 1,
                                 'stop': ref_stop + ref_start_phase})
        prefix = str(prefix).lower()

        # begin constructing reference and alternate codon sequences
        ref_codon = ref + prefix
        alt_codon = alt + prefix

        # obtain any subsequence nucleotides to complete the last codon
        if len(ref_codon) % 3:
            ref_stop_phase = len(ref_codon) % 3
            suffix = fasta.sequence({'chr': chrom,
                                     'start': ref_start - 3 + ref_stop_phase,
                                     'stop': ref_start - 1})
            suffix = str(suffix).lower()
            ref_codon = suffix + ref_codon

        if len(alt_codon) % 3:
            alt_stop_phase = len(alt_codon) % 3
            suffix = fasta.sequence({'chr': chrom,
                                     'start': ref_start - 3 + alt_stop_phase,
                                     'stop': ref_start - 1})
            suffix = str(suffix).lower()
            alt_codon = suffix + alt_codon

        # take reverse complement
        ref_codon = str(Seq(ref_codon).reverse_complement())
        alt_codon = str(Seq(alt_codon).reverse_complement())

    return ref_cds_start, ref_cds_stop, ref_start_phase, ref_codon, alt_codon


def _get_coding_position(ref_start, ref_stop, cds, cdss):

    if cds.strand == '+':

        # sort exons
        cdss = sorted(cdss, key=operator.itemgetter('start'))

        # find index of overlapping exons in all exons
        cds_ids = [f.feature_id for f in cdss]
        cds_index = cds_ids.index(cds.feature_id)

        # find offset
        offset = sum([f.stop - f.start + 1 for f in cdss[:cds_index]])

        # find ref cds position
        ref_cds_start = offset + (ref_start - cds.start)
        ref_cds_stop = offset + (ref_stop - cds.start)

    else:

        # sort exons (backwards this time)
        cdss = sorted(cdss, key=operator.itemgetter('stop'), reverse=True)

        # find index of overlapping exons in all exons
        cds_ids = [f.feature_id for f in cdss]
        cds_index = cds_ids.index(cds.feature_id)

        # find offset
        offset = sum([cds.stop - cds.start + 1 for cds in cdss[:cds_index]])

        # find ref cds position
        ref_cds_start = offset + (cds.stop - ref_stop)
        ref_cds_stop = offset + (cds.stop - ref_start)

    return ref_cds_start, ref_cds_stop

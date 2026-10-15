from geenuff.base.helpers import make_enum

# data-version of so.obo; terms of the older release 2.5.3 that SO renamed or made obsolete stay, as
# older GFF3 files still use them
_so_version = '2026-08-07'

HISTIDYL_TRNA = "histidyl_tRNA"
C_D_BOX_SNORNA_PRIMARY_TRANSCRIPT = "C_D_box_snoRNA_primary_transcript"
POLYCISTRONIC_MRNA = "polycistronic_mRNA"
RRNA_23S = "rRNA_23S"
TYROSINE_TRNA_PRIMARY_TRANSCRIPT = "tyrosine_tRNA_primary_transcript"
SHRNA_PRIMARY_TRANSCRIPT = "shRNA_primary_transcript"
CONTIG = "contig"
U3_SNORNA = "U3_snoRNA"
SUPERCONTIG = "supercontig"
CIRCULAR_SINGLE_STRANDED_DNA_CHROMOSOME = "circular_single_stranded_DNA_chromosome"
RNA_APTAMER = "RNA_aptamer"
ASX_TURN_RIGHT_HANDED_TYPE_ONE = "asx_turn_right_handed_type_one"
RECOMBINATION_HOTSPOT = "recombination_hotspot"
GENE_SEGMENT = "gene_segment"
J_GENE_RECOMBINATION_FEATURE = "J_gene_recombination_feature"
POLYPEPTIDE_LIGAND_CONTACT = "polypeptide_ligand_contact"
OXYS_RNA = "OxyS_RNA"
MRNA_WITH_PLUS_2_FRAMESHIFT = "mRNA_with_plus_2_frameshift"
LOW_COMPLEXITY_REGION = "low_complexity_region"
TARGETING_VECTOR = "targeting_vector"
GENETIC_MARKER = "genetic_marker"
CONSENSUS_AFLP_FRAGMENT = "consensus_AFLP_fragment"
BREU_MOTIF = "BREu_motif"
GENE_WITH_TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT = "gene_with_transcript_with_translational_frameshift"
COPY_NUMBER_VARIATION = "copy_number_variation"
EDITED_TRANSCRIPT = "edited_transcript"
INTERNAL_ELIMINATED_SEQUENCE = "internal_eliminated_sequence"
ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT = "engineered_foreign_transposable_element"
LYSOSOMAL_LOCALIZATION_SIGNAL = "lysosomal_localization_signal"
FIVE_PRIME_UST = "five_prime_UST"
BACTERIAL_RNAPOL_PROMOTER_SIGMA_70 = "bacterial_RNApol_promoter_sigma_70"
VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_REARRANGED_SEGMENT = "vertebrate_immunoglobulin_T_cell_receptor_rearranged_segment"
DSR_MOTIF = "DSR_motif"
FIVE_PRIME_INTRON = "five_prime_intron"
GENE_TRAP_CONSTRUCT = "gene_trap_construct"
LOR = "loR"
ABERRANT_PROCESSED_TRANSCRIPT = "aberrant_processed_transcript"
DMV4_MOTIF = "DMv4_motif"
V_NONAMER = "V_nonamer"
STEROL_REGULATORY_ELEMENT = "sterol_regulatory_element"
FIVE_PRIME_D_HEPTAMER = "five_prime_D_heptamer"
PEROXYWYBUTOSINE = "peroxywybutosine"
PRIMARY_TRANSCRIPT_REGION = "primary_transcript_region"
LAMBDA_VECTOR = "lambda_vector"
TWO_PRIME_O_METHYLPSEUDOURIDINE = "two_prime_O_methylpseudouridine"
INTERNAL_UTR = "internal_UTR"
TILING_PATH_CLONE = "tiling_path_clone"
BACTERIAL_RNAPOL_PROMOTER = "bacterial_RNApol_promoter"
DS_OLIGO = "ds_oligo"
TARGET_SITE_DUPLICATION = "target_site_duplication"
DPE_MOTIF = "DPE_motif"
SHADOW_ENHANCER = "shadow_enhancer"
CRM = "CRM"
NEGATIVE_SENSE_SSRNA_VIRAL_SEQUENCE = "negative_sense_ssRNA_viral_sequence"
HOMOL_E_BOX = "homol_E_box"
D_DJ_CLUSTER = "D_DJ_cluster"
FIVE_PRIME_D_SPACER = "five_prime_D_spacer"
THREE_TWO_PRIME_O_DIMETHYLURIDINE = "three_two_prime_O_dimethyluridine"
EST_MATCH = "EST_match"
GLUTAMINE = "glutamine"
METHIONYL_TRNA = "methionyl_tRNA"
PHENYLALANINE = "phenylalanine"
TRANSMEMBRANE_POLYPEPTIDE_REGION = "transmembrane_polypeptide_region"
CHROMOSOMAL_REGULATORY_ELEMENT = "chromosomal_regulatory_element"
ARGININE = "arginine"
ENDONUCLEASE_SPLICED_INTRON = "endonuclease_spliced_intron"
TRANSPOSABLE_ELEMENT_INSERTION_SITE = "transposable_element_insertion_site"
MACRONUCLEUS_DESTINED_SEGMENT = "macronucleus_destined_segment"
V_VJ_C_CLUSTER = "V_VJ_C_cluster"
MODIFIED_CYTOSINE = "modified_cytosine"
CODING_START = "coding_start"
TWO_METHYLTHIO_N6_HYDROXYNORVALYL_CARBAMOYLADENOSINE = "two_methylthio_N6_hydroxynorvalyl_carbamoyladenosine"
U3_THREE_PRIME_LTR_REGION = "U3_three_prime_LTR_region"
VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_GENE_CLUSTER = "vertebrate_immunoglobulin_T_cell_receptor_gene_cluster"
CIS_ACTING_HOMOLOGOUS_CHROMOSOME_PAIRING_REGION = "cis_acting_homologous_chromosome_pairing_region"
U3_LTR_REGION = "U3_LTR_region"
READING_FRAME = "reading_frame"
SUBTELOMERE = "subtelomere"
INVALIDATED_CDNA_CLONE = "invalidated_cDNA_clone"
SPLICED_LEADER_RNA = "spliced_leader_RNA"
TILING_PATH_FRAGMENT = "tiling_path_fragment"
IRLINV_SITE = "IRLinv_site"
PROVIRAL_REGION = "proviral_region"
TRANSPOSABLE_ELEMENT_FLANKING_REGION = "transposable_element_flanking_region"
TRANSCRIPTION_START_CLUSTER = "transcription_start_cluster"
V_GENE_SEGMENT = "V_gene_segment"
EXON_JUNCTION = "exon_junction"
SAP1_RECOGNITION_MOTIF = "Sap1_recognition_motif"
DISABLED_READING_FRAME = "disabled_reading_frame"
DEFECTIVE_CONJUGATIVE_TRANSPOSON = "defective_conjugative_transposon"
T_TO_A_TRANSVERSION = "T_to_A_transversion"
VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_SIGNAL_FEATURE = "vertebrate_immune_system_gene_recombination_signal_feature"
THREE_PRIME_UTR = "three_prime_UTR"
STRUCTURAL_ALTERATION = "structural_alteration"
CAPPED_MRNA = "capped_mRNA"
PROTEIN_CODING_GENE = "protein_coding_gene"
H3K36_MONOMETHYLATION_SITE = "H3K36_monomethylation_site"
_8_OXOADENINE = "8_oxoadenine"
SPLICING_REGULATORY_REGION = "splicing_regulatory_region"
HOMOL_D_BOX = "homol_D_box"
CENTROMERE_DNA_ELEMENT_I = "centromere_DNA_Element_I"
EST = "EST"
ENGINEERED_FOREIGN_GENE = "engineered_foreign_gene"
EXPRESSED_SEQUENCE_MATCH = "expressed_sequence_match"
TRANSCRIPT_BOUND_BY_PROTEIN = "transcript_bound_by_protein"
ATTB_SITE = "attB_site"
PRE_EDITED_MRNA = "pre_edited_mRNA"
DICISTRONIC_PRIMARY_TRANSCRIPT = "dicistronic_primary_transcript"
MODIFIED_L_METHIONINE = "modified_L_methionine"
UAA_STOP_CODON_SIGNAL = "UAA_stop_codon_signal"
NATURAL_TRANSPOSABLE_ELEMENT = "natural_transposable_element"
FIVE_PRIME_LTR_COMPONENT = "five_prime_LTR_component"
PAC = "PAC"
P_ELEMENT = "p_element"
CENTROMERE_DNA_ELEMENT_III = "centromere_DNA_Element_III"
LINCRNA = "lincRNA"
BETA_BULGE_LOOP = "beta_bulge_loop"
RPRA_RNA = "RprA_RNA"
ANCHOR_BINDING_SITE = "anchor_binding_site"
D_DJ_J_C_CLUSTER = "D_DJ_J_C_cluster"
TWO_PRIME_O_METHYLURIDINE = "two_prime_O_methyluridine"
Y_RNA_PRIMARY_TRANSCRIPT = "Y_RNA_primary_transcript"
A_MINOR_RNA_MOTIF = "A_minor_RNA_motif"
AU_RICH_ELEMENT = "AU_rich_element"
FIVE_PRIME_OPEN_READING_FRAME = "five_prime_open_reading_frame"
UNIGENE_CLUSTER = "unigene_cluster"
POLYPEPTIDE_REPEAT = "polypeptide_repeat"
RNA_JUNCTION_LOOP = "RNA_junction_loop"
COINTEGRATED_PLASMID = "cointegrated_plasmid"
CLASS_II_RNA = "class_II_RNA"
IRRINV_SITE = "IRRinv_site"
RRNA_LARGE_SUBUNIT_PRIMARY_TRANSCRIPT = "rRNA_large_subunit_primary_transcript"
CDNA_MATCH = "cDNA_match"
NDM3_MOTIF = "NDM3_motif"
CLONE_INSERT_END = "clone_insert_end"
CENTROMERE = "centromere"
SELENOCYSTEINYL_TRNA = "selenocysteinyl_tRNA"
DG_REPEAT = "dg_repeat"
PROLINE_TRNA_PRIMARY_TRANSCRIPT = "proline_tRNA_primary_transcript"
TRANSCRIBED_SPACER_REGION = "transcribed_spacer_region"
INTEGRATED_PLASMID = "integrated_plasmid"
CATMAT_RIGHT_HANDED_THREE = "catmat_right_handed_three"
REGIONAL_CENTROMERE_OUTER_REPEAT_TRANSCRIPT = "regional_centromere_outer_repeat_transcript"
SEVEN_METHYLGUANOSINE = "seven_methylguanosine"
C_TO_T_TRANSITION = "C_to_T_transition"
POLYPEPTIDE_DNA_CONTACT = "polypeptide_DNA_contact"
NON_TERMINAL_RESIDUE = "non_terminal_residue"
SEQUENCE_DIFFERENCE = "sequence_difference"
REPEAT_FRAGMENT = "repeat_fragment"
DMV3_MOTIF = "DMv3_motif"
TSS_REGION = "TSS_region"
RESTRICTION_ENZYME_SINGLE_STRAND_OVERHANG = "restriction_enzyme_single_strand_overhang"
ENHANCER_TRAP_CONSTRUCT = "enhancer_trap_construct"
PSEUDOGENE_BY_UNEQUAL_CROSSING_OVER = "pseudogene_by_unequal_crossing_over"
ALPHA_BETA_MOTIF = "alpha_beta_motif"
APICOPLAST_CHROMOSOME = "apicoplast_chromosome"
N6_ACETYLADENOSINE = "N6_acetyladenosine"
MODIFIED_L_GLUTAMIC_ACID = "modified_L_glutamic_acid"
FIVE_PRIME_CIS_SPLICE_SITE = "five_prime_cis_splice_site"
SYNTHETIC_OLIGO = "synthetic_oligo"
NON_ADJACENT_RESIDUES = "non_adjacent_residues"
CAP = "cap"
ASSEMBLY = "assembly"
CYSTEINE = "cysteine"
GROUP_II_INTRON = "group_II_intron"
RECOMBINATION_REGULATORY_REGION = "recombination_regulatory_region"
MODIFIED_L_ASPARAGINE = "modified_L_asparagine"
SECIS_ELEMENT = "SECIS_element"
U5_THREE_PRIME_LTR_REGION = "U5_three_prime_LTR_region"
REVERSE_HOOGSTEEN_BASE_PAIR = "reverse_Hoogsteen_base_pair"
OLIGO = "oligo"
SMFISH_PROBE = "smFISH_probe"
H_ACA_BOX_SNORNA = "H_ACA_box_snoRNA"
V_D_DJ_CLUSTER = "V_D_DJ_cluster"
PARTIAL_GENOMIC_SEQUENCE_ASSEMBLY = "partial_genomic_sequence_assembly"
PSEUDOGENIC_REGION = "pseudogenic_region"
CONSENSUS_MRNA = "consensus_mRNA"
HELITRON = "helitron"
SEVEN_DEAZAGUANOSINE = "seven_deazaguanosine"
THREE_PRIME_LTR = "three_prime_LTR"
CHROMOSOME_BREAKPOINT = "chromosome_breakpoint"
SNRNA = "snRNA"
SNORNA_GENE = "snoRNA_gene"
H4K_ACYLATION_REGION = "H4K_acylation_region"
I_BOX = "I-box"
H3K27_TRIMETHYLATION_SITE = "H3K27_trimethylation_site"
TRANS_SPLICE_ACCEPTOR_SITE = "trans_splice_acceptor_site"
LINEAR_SINGLE_STRANDED_RNA_CHROMOSOME = "linear_single_stranded_RNA_chromosome"
ASPE_PRIMER = "ASPE_primer"
PNA_OLIGO = "PNA_oligo"
MITE = "MITE"
SEQUENCE_MOTIF = "sequence_motif"
ASPARAGINYL_TRNA = "asparaginyl_tRNA"
THREONYL_TRNA = "threonyl_tRNA"
TYROSYL_TRNA = "tyrosyl_tRNA"
SELENOCYSTEINE = "selenocysteine"
MRNA_WITH_MINUS_2_FRAMESHIFT = "mRNA_with_minus_2_frameshift"
GENE_GROUP_REGULATORY_REGION = "gene_group_regulatory_region"
ENHANCER_BINDING_SITE = "enhancer_binding_site"
MRNA_RECODED_BY_TRANSLATIONAL_BYPASS = "mRNA_recoded_by_translational_bypass"
BASE_CALL_ERROR_CORRECTION = "base_call_error_correction"
CDS_INDEPENDENTLY_KNOWN = "CDS_independently_known"
BETA_TURN_LEFT_HANDED_TYPE_TWO = "beta_turn_left_handed_type_two"
GAMMA_TURN = "gamma_turn"
RESTRICTION_FRAGMENT = "restriction_fragment"
SP6_RNA_POLYMERASE_PROMOTER = "SP6_RNA_Polymerase_Promoter"
TRANSCRIPT_REGION = "transcript_region"
SYMMETRIC_RNA_INTERNAL_LOOP = "symmetric_RNA_internal_loop"
MINUS_24_SIGNAL = "minus_24_signal"
ONE_TWO_PRIME_O_DIMETHYLINOSINE = "one_two_prime_O_dimethylinosine"
POLYPEPTIDE_MANGANESE_ION_CONTACT_SITE = "polypeptide_manganese_ion_contact_site"
CHIMERIC_CDNA_CLONE = "chimeric_cDNA_clone"
PARTIALLY_PROCESSED_CDNA_CLONE = "partially_processed_cDNA_clone"
INVERSION = "inversion"
HEPTAMER_OF_RECOMBINATION_FEATURE_OF_VERTEBRATE_IMMUNE_SYSTEM_GENE = "heptamer_of_recombination_feature_of_vertebrate_immune_system_gene"
FIVE_HYDROXYURIDINE = "five_hydroxyuridine"
V_J_C_CLUSTER = "V_J_C_cluster"
H3K9_DIMETHYLATION_SITE = "H3K9_dimethylation_site"
CSL_RESPONSE_ELEMENT = "CSL_response_element"
ST_TURN_LEFT_HANDED_TYPE_TWO = "st_turn_left_handed_type_two"
CRE = "CRE"
POLYPEPTIDE_MOLYBDENUM_ION_CONTACT_SITE = "polypeptide_molybdenum_ion_contact_site"
G_BOX = "G_box"
TRANS_SPLICE_JUNCTION = "trans_splice_junction"
HELIX_TURN_HELIX = "helix_turn_helix"
GENE_GROUP = "gene_group"
PRIMER_MATCH = "primer_match"
PRIRNA = "priRNA"
FIVE_PRIME_UTR_INTRON = "five_prime_UTR_intron"
CDS_SUPPORTED_BY_SEQUENCE_SIMILARITY_DATA = "CDS_supported_by_sequence_similarity_data"
CHROMOSOME = "chromosome"
NUCLEOMORPH_GENE = "nucleomorph_gene"
_5_METHYLCYTOSINE = "5_methylcytosine"
SL3_ACCEPTOR_SITE = "SL3_acceptor_site"
H4K4_TRIMETHYLATION_SITE = "H4K4_trimethylation_site"
INTRONIC_SPLICING_SILENCER = "intronic_splicing_silencer"
N6_2_PRIME_O_DIMETHYLADENOSINE = "N6_2_prime_O_dimethyladenosine"
METAL_BINDING_SITE = "metal_binding_site"
DMV5_MOTIF = "DMv5_motif"
POLYPEPTIDE_BINDING_MOTIF = "polypeptide_binding_motif"
GTT_TRINUCLEOTIDE_REPEAT = "GTT_trinucleotide_repeat"
U5_LTR_REGION = "U5_LTR_region"
N6_METHYLADENINE = "N6_methyladenine"
PYRIMIDINE_TRANSITION = "pyrimidine_transition"
DH_REPEAT = "dh_repeat"
GLUTAMIC_ACID_TRNA_PRIMARY_TRANSCRIPT = "glutamic_acid_tRNA_primary_transcript"
CONTIG_READ = "contig_read"
PRE_EDITED_REGION = "pre_edited_region"
TRANSPOSABLE_ELEMENT_CDS = "transposable_element_CDS"
H4R3_DIMETHYLATION_SITE = "H4R3_dimethylation_site"
CATMAT_LEFT_HANDED_FOUR = "catmat_left_handed_four"
ISOWYOSINE = "isowyosine"
SS_RNA_VIRAL_SEQUENCE = "ss_RNA_viral_sequence"
N6_N6_DIMETHYLADENOSINE = "N6_N6_dimethyladenosine"
VACUOLAR_SORTING_SIGNAL = "vacuolar_sorting_signal"
PROTEIN_CODING_PRIMARY_TRANSCRIPT = "protein_coding_primary_transcript"
SNRNA_PRIMARY_TRANSCRIPT = "snRNA_primary_transcript"
TRANSCRIPTION_REGULATORY_REGION = "transcription_regulatory_region"
MIRNA_LOOP = "miRNA_loop"
DOUBLE_STRANDED_DNA_CHROMOSOME = "double_stranded_DNA_chromosome"
INR1_MOTIF = "INR1_motif"
L_BOX = "L_box"
N6_THREONYLCARBAMOYLADENOSINE = "N6_threonylcarbamoyladenosine"
SEVEN_AMINOMETHYL_SEVEN_DEAZAGUANOSINE = "seven_aminomethyl_seven_deazaguanosine"
ONE_METHYLGUANOSINE = "one_methylguanosine"
FIVE_PRIME_UTR = "five_prime_UTR"
UST_MATCH = "UST_match"
DCE_SII = "DCE_SII"
RR_TRACT = "RR_tract"
EPIGENETICALLY_MODIFIED_REGION = "epigenetically_modified_region"
TELOMERASE_RNA_GENE = "telomerase_RNA_gene"
PATERNALLY_IMPRINTED_GENE = "paternally_imprinted_gene"
PSEUDOGENIC_TRNA = "pseudogenic_tRNA"
PROPHAGE = "prophage"
RNA_6S = "RNA_6S"
V_VDJ_C_CLUSTER = "V_VDJ_C_cluster"
INTERBAND = "interband"
INSERTION = "insertion"
FIVE_ISOPENTENYLAMINOMETHYL_TWO_PRIME_O_METHYLURIDINE = "five_isopentenylaminomethyl_two_prime_O_methyluridine"
THREE_METHYLURIDINE = "three_methyluridine"
H3K27_ACETYLATION_SITE = "H3K27_acetylation_site"
CIRCULAR_DOUBLE_STRANDED_DNA_CHROMOSOME = "circular_double_stranded_DNA_chromosome"
LINCRNA_GENE = "lincRNA_gene"
CHROMOSOME_BREAKAGE_SEQUENCE = "chromosome_breakage_sequence"
CTG_START_CODON = "CTG_start_codon"
H3K18_ACETYLATION_SITE = "H3K18_acetylation_site"
POLYPEPTIDE = "polypeptide"
H4K8_ACYLATION_SITE = "H4K8_acylation_site"
STOP_CODON_SIGNAL = "stop_codon_signal"
DELETION = "deletion"
THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = "three_prime_sticky_end_restriction_enzyme_cleavage_site"
SIGNATURE = "signature"
LEUCINE_TRNA_PRIMARY_TRANSCRIPT = "leucine_tRNA_primary_transcript"
MT_GENE = "mt_gene"
RESCUE_REGION = "rescue_region"
NESTED_REPEAT = "nested_repeat"
MONOCISTRONIC_MRNA = "monocistronic_mRNA"
H3K36_DIMETHYLATION_SITE = "H3K36_dimethylation_site"
GENE = "gene"
MULTIPLEXING_SEQUENCE_IDENTIFIER = "multiplexing_sequence_identifier"
MODIFIED_L_HISTIDINE = "modified_L_histidine"
X_ELEMENT = "X_element"
BACTERIAL_RNAPOL_PROMOTER_SIGMA_ECF = "bacterial_RNApol_promoter_sigma_ecf"
TRANSLATIONALLY_REGULATED_GENE = "translationally_regulated_gene"
METHYLWYOSINE = "methylwyosine"
ENGINEERED_PLASMID = "engineered_plasmid"
ENGINEERED_EPISOME = "engineered_episome"
SINGLE_STRAND_RESTRICTION_ENZYME_CLEAVAGE_SITE = "single_strand_restriction_enzyme_cleavage_site"
CDNA_CLONE = "cDNA_clone"
MODIFIED_L_PHENYLALANINE = "modified_L_phenylalanine"
CHROMOSOME_ARM = "chromosome_arm"
TWO_METHYLTHIO_N6_THREONYL_CARBAMOYLADENOSINE = "two_methylthio_N6_threonyl_carbamoyladenosine"
NON_TRANSCRIBED_REGION = "non_transcribed_region"
UNEDITED_REGION = "unedited_region"
BLUNT_END_RESTRICTION_ENZYME_CLEAVAGE_JUNCTION = "blunt_end_restriction_enzyme_cleavage_junction"
EXON_REGION = "exon_region"
RIBOSOME_ENTRY_SITE = "ribosome_entry_site"
POLYPYRIMIDINE_TRACT = "polypyrimidine_tract"
SEQUENCE_LENGTH_VARIATION = "sequence_length_variation"
THREE_PRIME_RACE_CLONE = "three_prime_RACE_clone"
GENOMIC_CLONE = "genomic_clone"
EDITED_CDS = "edited_CDS"
POLYA_SEQUENCE = "polyA_sequence"
UNDERMODIFIED_HYDROXYWYBUTOSINE = "undermodified_hydroxywybutosine"
FRE = "FRE"
ASX_TURN = "asx_turn"
MIR_ENCODING_LNCRNA_PRIMARY_TRANSCRIPT = "miR_encoding_lncRNA_primary_transcript"
THREE_PRIME_RST = "three_prime_RST"
THREE_PRIME_D_SPACER = "three_prime_D_spacer"
CATALYTIC_RESIDUE = "catalytic_residue"
SILENCER = "silencer"
SL6_ACCEPTOR_SITE = "SL6_acceptor_site"
AMBISENSE_SSRNA_VIRAL_SEQUENCE = "ambisense_ssRNA_viral_sequence"
TRANS_SPLICE_DONOR_SITE = "trans_splice_donor_site"
TMRNA_REGION = "tmRNA_region"
BETA_TURN_TYPE_SIX_A_TWO = "beta_turn_type_six_a_two"
ONE_TWO_PRIME_O_DIMETHYLGUANOSINE = "one_two_prime_O_dimethylguanosine"
PROTEIN_STABILITY_ELEMENT = "protein_stability_element"
NUCLEAR_LOCALIZATION_SIGNAL = "nuclear_localization_signal"
MATURE_TRANSCRIPT_REGION = "mature_transcript_region"
MRNA_WITH_PLUS_1_FRAMESHIFT = "mRNA_with_plus_1_frameshift"
ARS_CONSENSUS_SEQUENCE = "ARS_consensus_sequence"
CIS_SPLICE_SITE = "cis_splice_site"
EXPRESSED_SEQUENCE_ASSEMBLY = "expressed_sequence_assembly"
EXON = "exon"
DISTAL_DUPLICATION = "distal_duplication"
MIRNA_STEM = "miRNA_stem"
LNCRNA_PRIMARY_TRANSCRIPT = "lncRNA_primary_transcript"
ENZYMATIC_RNA = "enzymatic_RNA"
TWO_PRIME_O_RIBOSYLADENOSINE_PHOSPHATE = "two_prime_O_ribosyladenosine_phosphate"
TRNA_GENE = "tRNA_gene"
CAAT_SIGNAL = "CAAT_signal"
POLYPEPTIDE_STRUCTURAL_REGION = "polypeptide_structural_region"
ATTENUATOR = "attenuator"
C_TERMINAL_REGION = "c_terminal_region"
ORPHAN_CDS = "orphan_CDS"
GRNA_GENE = "gRNA_gene"
LINKAGE_GROUP = "linkage_group"
SEQUENCE_UNCERTAINTY = "sequence_uncertainty"
U14_SNORNA_PRIMARY_TRANSCRIPT = "U14_snoRNA_primary_transcript"
SAGE_TAG = "SAGE_tag"
TRANSLOCATION = "translocation"
COMPLEX_SUBSTITUTION = "complex_substitution"
V_DJ_J_C_CLUSTER = "V_DJ_J_C_cluster"
FIVE_ISOPENTENYLAMINOMETHYL_URIDINE = "five_isopentenylaminomethyl_uridine"
LTR_COMPONENT = "LTR_component"
DJ_GENE_SEGMENT = "DJ_gene_segment"
THREONINE_TRNA_PRIMARY_TRANSCRIPT = "threonine_tRNA_primary_transcript"
RST = "RST"
N6_METHYLADENOSINE = "N6_methyladenosine"
N6_CIS_HYDROXYISOPENTENYL_ADENOSINE = "N6_cis_hydroxyisopentenyl_adenosine"
FIVE_HYDROXYMETHYLCYTIDINE = "five_hydroxymethylcytidine"
NUCLEAR_GENE = "nuclear_gene"
V_HEPTAMER = "V_heptamer"
BETA_TURN_RIGHT_HANDED_TYPE_TWO = "beta_turn_right_handed_type_two"
MODIFIED_DNA_BASE = "modified_DNA_base"
INR_MOTIF = "INR_motif"
N4_ACETYLCYTIDINE = "N4_acetylcytidine"
PROPLASTID_GENE = "proplastid_gene"
SEQUENCING_PRIMER = "sequencing_primer"
ONE_METHYLPSEUDOURIDINE = "one_methylpseudouridine"
MODIFIED_AMINO_ACID_FEATURE = "modified_amino_acid_feature"
BETA_BULGE_LOOP_FIVE = "beta_bulge_loop_five"
HAMMERHEAD_RIBOZYME = "hammerhead_ribozyme"
MCB = "MCB"
NCRNA_GENE = "ncRNA_gene"
Z1_REGION = "Z1_region"
H_PSEUDOKNOT = "H_pseudoknot"
OPEN_CHROMATIN_REGION = "open_chromatin_region"
MIR_ENCODING_VAULTRNA_PRIMARY_TRANSCRIPT = "miR_encoding_vaultRNA_primary_transcript"
TRANSPOSABLE_ELEMENT_PSEUDOGENE = "transposable_element_pseudogene"
ZINC_REPRESSED_ELEMENT = "zinc_repressed_element"
POSITIVELY_AUTOREGULATED_GENE = "positively_autoregulated_gene"
D_J_C_CLUSTER = "D_J_C_cluster"
LEUCOPLAST_GENE = "leucoplast_gene"
TMRNA_GENE = "tmRNA_gene"
INTRACHROMOSOMAL_BREAKPOINT = "intrachromosomal_breakpoint"
ENGINEERED_TAG = "engineered_tag"
WYOSINE = "wyosine"
GENE_CASSETTE = "gene_cassette"
U2_SNRNA = "U2_snRNA"
NON_LTR_RETROTRANSPOSON_POLYMERIC_TRACT = "non_LTR_retrotransposon_polymeric_tract"
POLYCISTRONIC_TRANSCRIPT = "polycistronic_transcript"
BETA_TURN_TYPE_EIGHT = "beta_turn_type_eight"
MODIFIED_L_GLUTAMINE = "modified_L_glutamine"
N6_GLYCINYLCARBAMOYLADENOSINE = "N6_glycinylcarbamoyladenosine"
H2AK9_ACETYLATION_SITE = "H2AK9_acetylation_site"
LNA_OLIGO = "LNA_oligo"
N_REGION = "N_region"
CAGE_TAG = "CAGE_tag"
THREE_PRIME_EST = "three_prime_EST"
GENE_WITH_RECODED_MRNA = "gene_with_recoded_mRNA"
BIOCHEMICAL_REGION_OF_PEPTIDE = "biochemical_region_of_peptide"
INTEGRATED_MOBILE_GENETIC_ELEMENT = "integrated_mobile_genetic_element"
Y_RNA = "Y_RNA"
TRANS_SPLICED_TRANSCRIPT = "trans_spliced_transcript"
WC_BASE_PAIR = "WC_base_pair"
ARCHAEOSINE = "archaeosine"
PROMOTER = "promoter"
ORIC = "oriC"
GENE_REARRANGED_AT_DNA_LEVEL = "gene_rearranged_at_DNA_level"
DDB_BOX = "DDB_box"
H2BK120_ACETYLATION_SITE = "H2BK120_acetylation_site"
LEUCYL_TRNA = "leucyl_tRNA"
TRNA_INTRON = "tRNA_intron"
OPERON = "operon"
H3K4AC_ACETYLATION_SITE = "H3K4ac_acetylation_site"
ER_RETENTION_SIGNAL = "ER_retention_signal"
FIVE_TWO_PRIME_O_DIMETHYLCYTIDINE = "five_two_prime_O_dimethylcytidine"
TMRNA_ACCEPTOR_PIECE = "tmRNA_acceptor_piece"
PIRNA_GENE = "piRNA_gene"
POLYPEPTIDE_COPPER_ION_CONTACT_SITE = "polypeptide_copper_ion_contact_site"
GLUTAMYL_TRNA = "glutamyl_tRNA"
RNA_STABILITY_ELEMENT = "RNA_stability_element"
NUCLEASE_SENSITIVE_SITE = "nuclease_sensitive_site"
CLIP = "clip"
_4_METHYLCYTOSINE = "4_methylcytosine"
ARGININE_TRNA_PRIMARY_TRANSCRIPT = "arginine_tRNA_primary_transcript"
FLANKING_THREE_PRIME_QUADRUPLET_RECODING_SIGNAL = "flanking_three_prime_quadruplet_recoding_signal"
MATURE_TRANSCRIPT = "mature_transcript"
TRANSLATED_NUCLEOTIDE_MATCH = "translated_nucleotide_match"
SUBSTITUTION = "substitution"
GENE_FRAGMENT = "gene_fragment"
GENE_WITH_STOP_CODON_REDEFINED_AS_PYRROLYSINE = "gene_with_stop_codon_redefined_as_pyrrolysine"
MITOCHONDRIAL_CHROMOSOME = "mitochondrial_chromosome"
MODIFIED_L_ALANINE = "modified_L_alanine"
REGIONAL_CENTROMERE_CENTRAL_CORE = "regional_centromere_central_core"
MATING_TYPE_REGION = "mating_type_region"
NUCLEOTIDE_TO_PROTEIN_BINDING_SITE = "nucleotide_to_protein_binding_site"
MODIFIED_URIDINE = "modified_uridine"
RFLP_FRAGMENT = "RFLP_fragment"
CYTOPLASMIC_POLYPEPTIDE_REGION = "cytoplasmic_polypeptide_region"
D_DJ_C_CLUSTER = "D_DJ_C_cluster"
TRANSMEMBRANE_HELIX = "transmembrane_helix"
V_REGION = "V_region"
GAP = "gap"
V_J_CLUSTER = "V_J_cluster"
PCB = "PCB"
INTEGRATION_EXCISION_SITE = "integration_excision_site"
LOCUS_CONTROL_REGION = "locus_control_region"
C_REGION = "C_region"
DNA_CONSTRAINT_SEQUENCE = "DNA_constraint_sequence"
LONG_TERMINAL_REPEAT = "long_terminal_repeat"
MATRIX_ATTACHMENT_SITE = "matrix_attachment_site"
TWO_THIOCYTIDINE = "two_thiocytidine"
REVERSE_PRIMER = "reverse_primer"
GLYCINE = "glycine"
TELOMERIC_REPEAT = "telomeric_repeat"
THREE_PRIME_CODING_EXON_CODING_REGION = "three_prime_coding_exon_coding_region"
MINICIRCLE_GENE = "minicircle_gene"
POLYPEPTIDE_NEST_MOTIF = "polypeptide_nest_motif"
NON_CANONICAL_THREE_PRIME_SPLICE_SITE = "non_canonical_three_prime_splice_site"
CROSS_GENOME_MATCH = "cross_genome_match"
ACE2_UAS = "Ace2_UAS"
CDRE_MOTIF = "CDRE_motif"
FIVE_ISOPENTENYLAMINOMETHYL_TWO_THIOURIDINE = "five_isopentenylaminomethyl_two_thiouridine"
POST_TRANSLATIONALLY_REGULATED_GENE = "post_translationally_regulated_gene"
PLASMID_VECTOR = "plasmid_vector"
VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_SPACER = "vertebrate_immune_system_gene_recombination_spacer"
TRANSCRIBED_FRAGMENT = "transcribed_fragment"
GENE_WITH_EDITED_TRANSCRIPT = "gene_with_edited_transcript"
RECOMBINATIONALLY_REARRANGED_GENE = "recombinationally_rearranged_gene"
ORIGIN_OF_REPLICATION = "origin_of_replication"
DMV1_MOTIF = "DMv1_motif"
MINUS_10_SIGNAL = "minus_10_signal"
VJ_GENE_SEGMENT = "VJ_gene_segment"
PROCESSED_TRANSCRIPT = "processed_transcript"
NUCLEASE_HYPERSENSITIVE_SITE = "nuclease_hypersensitive_site"
TCT_MOTIF = "TCT_motif"
SEVEN_CYANO_SEVEN_DEAZAGUANOSINE = "seven_cyano_seven_deazaguanosine"
TF_BINDING_SITE = "TF_binding_site"
MNV = "MNV"
V_VJ_J_C_CLUSTER = "V_VJ_J_C_cluster"
RECOMBINATION_FEATURE_OF_REARRANGED_GENE = "recombination_feature_of_rearranged_gene"
GATA_BOX = "GATA_box"
RNAI_REAGENT = "RNAi_reagent"
REPEAT_UNIT = "repeat_unit"
SORF = "sORF"
TNA_OLIGO = "TNA_oligo"
PLASMID = "plasmid"
H3K27_MONOMETHYLATION_SITE = "H3K27_monomethylation_site"
REGIONAL_CENTROMERE_OUTER_REPEAT_REGION = "regional_centromere_outer_repeat_region"
MAJOR_TSS = "major_TSS"
MODIFIED_ADENOSINE = "modified_adenosine"
GENE_WITH_STOP_CODON_REDEFINED_AS_SELENOCYSTEINE = "gene_with_stop_codon_redefined_as_selenocysteine"
FORWARD_PRIMER = "forward_primer"
UNASSIGNED_SUPERCONTIG = "unassigned_supercontig"
GUIDE_RNA_REGION = "guide_RNA_region"
READ = "read"
VDJ_C_CLUSTER = "VDJ_C_cluster"
IMMATURE_PEPTIDE_REGION = "immature_peptide_region"
TATA_BOX = "TATA_box"
BAC = "BAC"
ASPARAGINE = "asparagine"
LIGATION_BASED_READ = "ligation_based_read"
NONAMER_OF_RECOMBINATION_FEATURE_OF_VERTEBRATE_IMMUNE_SYSTEM_GENE = "nonamer_of_recombination_feature_of_vertebrate_immune_system_gene"
H4K20_MONOMETHYLATION_SITE = "H4K20_monomethylation_site"
PAC_END = "PAC_end"
MANNOSYL_QUEUOSINE = "mannosyl_queuosine"
UORF = "uORF"
_5_CARBOXYLCYTOSINE = "5_carboxylcytosine"
POLYPEPTIDE_METAL_CONTACT = "polypeptide_metal_contact"
CDS_FRAGMENT = "CDS_fragment"
BETA_BULGE_LOOP_SIX = "beta_bulge_loop_six"
MIRNA_TARGET_SITE = "miRNA_target_site"
TRANS_SPLICE_SITE = "trans_splice_site"
RNA_MOTIF = "RNA_motif"
MODIFIED_L_THREONINE = "modified_L_threonine"
V_D_J_CLUSTER = "V_D_J_cluster"
SATELLITE_DNA = "satellite_DNA"
MODIFIED_INOSINE = "modified_inosine"
PSEUDOURIDINE = "pseudouridine"
VDJ_J_CLUSTER = "VDJ_J_cluster"
GAGA_MOTIF = "GAGA_motif"
PRIMER_BINDING_SITE = "primer_binding_site"
N6_METHYL_N6_THREONYLCARBAMOYLADENOSINE = "N6_methyl_N6_threonylcarbamoyladenosine"
RRNA_28S = "rRNA_28S"
DNA_MOTIF = "DNA_motif"
MICROARRAY_OLIGO = "microarray_oligo"
TRANSVERSION = "transversion"
DNA_CHROMOSOME = "DNA_chromosome"
INSULATOR_BINDING_SITE = "insulator_binding_site"
RECODING_STIMULATORY_REGION = "recoding_stimulatory_region"
POSSIBLE_ASSEMBLY_ERROR = "possible_assembly_error"
DISTANT_THREE_PRIME_RECODING_SIGNAL = "distant_three_prime_recoding_signal"
PROVIRAL_GENE = "proviral_gene"
PHAGEMID = "phagemid"
MACRONUCLEAR_CHROMOSOME = "macronuclear_chromosome"
V_D_DJ_J_C_CLUSTER = "V_D_DJ_J_C_cluster"
MOR = "moR"
TMRNA_PRIMARY_TRANSCRIPT = "tmRNA_primary_transcript"
THREE_PRIME_CIS_SPLICE_SITE = "three_prime_cis_splice_site"
U6_SNRNA = "U6_snRNA"
SEQUENCE_FEATURE = "sequence_feature"
DNAZYME = "DNAzyme"
ALPHA_HELIX = "alpha_helix"
HYDROPHOBIC_REGION_OF_PEPTIDE = "hydrophobic_region_of_peptide"
RNA_HOOK_TURN = "RNA_hook_turn"
LNC_RNA = "lncRNA"
STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = "sticky_end_restriction_enzyme_cleavage_site"
GENE_MEMBER_REGION = "gene_member_region"
SERINE_TRNA_PRIMARY_TRANSCRIPT = "serine_tRNA_primary_transcript"
CORE_PROMOTER_ELEMENT = "core_promoter_element"
PSEUDOGENIC_RRNA = "pseudogenic_rRNA"
RNA_POLYMERASE_III_TATA_BOX = "RNA_polymerase_III_TATA_box"
DJ_J_CLUSTER = "DJ_J_cluster"
POLYPEPTIDE_TUNGSTEN_ION_CONTACT_SITE = "polypeptide_tungsten_ion_contact_site"
SS_OLIGO = "ss_oligo"
U11_SNRNA = "U11_snRNA"
VJ_C_CLUSTER = "VJ_C_cluster"
ACTIVE_PEPTIDE = "active_peptide"
CYSTEINYL_TRNA = "cysteinyl_tRNA"
ST_TURN_RIGHT_HANDED_TYPE_ONE = "st_turn_right_handed_type_one"
CONFORMATIONAL_SWITCH = "conformational_switch"
SELENOCYSTEINE_TRNA_PRIMARY_TRANSCRIPT = "selenocysteine_tRNA_primary_transcript"
PEPTIDE_LOCALIZATION_SIGNAL = "peptide_localization_signal"
THREE_PRIME_LTR_COMPONENT = "three_prime_LTR_component"
BIOMATERIAL_REGION = "biomaterial_region"
GENE_WITH_NON_CANONICAL_START_CODON = "gene_with_non_canonical_start_codon"
READ_PAIR = "read_pair"
DINUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE = "dinucleotide_repeat_microsatellite_feature"
INVERTED_TANDEM_DUPLICATION = "inverted_tandem_duplication"
GOLDEN_PATH = "golden_path"
TRANSIT_PEPTIDE = "transit_peptide"
NOVEL_SEQUENCE_INSERTION = "novel_sequence_insertion"
RRNA_PRIMARY_TRANSCRIPT = "rRNA_primary_transcript"
CDS = "CDS"
MATURE_PROTEIN_REGION = "mature_protein_region"
EPIGENETICALLY_MODIFIED_GENE = "epigenetically_modified_gene"
H3K9_METHYLATION_SITE = "H3K9_methylation_site"
FOUR_BP_START_CODON = "four_bp_start_codon"
GENE_WITH_START_CODON_CUG = "gene_with_start_codon_CUG"
MUTATED_VARIANT_SITE = "mutated_variant_site"
MINUS_35_SIGNAL = "minus_35_signal"
ARS = "ARS"
ASSEMBLY_COMPONENT = "assembly_component"
TRANSPOSABLE_ELEMENT = "transposable_element"
J_HEPTAMER = "J_heptamer"
CT_GENE = "ct_gene"
NUCLEAR_CHROMOSOME = "nuclear_chromosome"
H2BK15_ACETYLATION_SITE = "H2BK15_acetylation_site"
SERINE_THREONINE_STAPLE_MOTIF = "serine_threonine_staple_motif"
STS_MAP = "STS_map"
CYANELLE_CHROMOSOME = "cyanelle_chromosome"
TWINTRON = "twintron"
CRYPTOGENE = "cryptogene"
KINETOPLAST_GENE = "kinetoplast_gene"
DS_DNA_VIRAL_SEQUENCE = "ds_DNA_viral_sequence"
HISTONE_BINDING_SITE = "histone_binding_site"
BACTERIAL_TERMINATOR = "bacterial_terminator"
POLYMERASE_SYNTHESIS_READ = "polymerase_synthesis_read"
TWO_METHYLTHIO_N6_METHYLADENOSINE = "two_methylthio_N6_methyladenosine"
U2_INTRON = "U2_intron"
FUNCTIONAL_CANDIDATE_GENE = "functional_candidate_gene"
HISTONE_UBIQITINATION_SITE = "histone_ubiqitination_site"
THREE_PRIME_D_NONAMER = "three_prime_D_nonamer"
NCRNA = "ncRNA"
H2BK12_ACETYLATION_SITE = "H2BK12_acetylation_site"
FIVE_CARBOXYHYDROXYMETHYL_URIDINE = "five_carboxyhydroxymethyl_uridine"
RETROTRANSPOSON = "retrotransposon"
INTRONIC_SPLICE_ENHANCER = "intronic_splice_enhancer"
SPECIFIC_RECOMBINATION_SITE = "specific_recombination_site"
CANDIDATE_GENE = "candidate_gene"
RRE_RNA = "RRE_RNA"
M26_BINDING_SITE = "M26_binding_site"
ISOLEUCYL_TRNA = "isoleucyl_tRNA"
A_TO_T_TRANSVERSION = "A_to_T_transversion"
HISTIDINE_TRNA_PRIMARY_TRANSCRIPT = "histidine_tRNA_primary_transcript"
MINISATELLITE = "minisatellite"
FIVE_CARBOXYMETHYLAMINOMETHYL_TWO_PRIME_O_METHYLURIDINE = "five_carboxymethylaminomethyl_two_prime_O_methyluridine"
SPOT_42_RNA = "spot_42_RNA"
H3K79_DIMETHYLATION_SITE = "H3K79_dimethylation_site"
DISTAL_PROMOTER_ELEMENT = "distal_promoter_element"
MODIFIED_L_SERINE = "modified_L_serine"
POLYPEPTIDE_CALCIUM_ION_CONTACT_SITE = "polypeptide_calcium_ion_contact_site"
CCAAT_MOTIF = "CCAAT_motif"
METABOLIC_ISLAND = "metabolic_island"
VAULTRNA_PRIMARY_TRANSCRIPT = "vaultRNA_primary_transcript"
RNAPOL_III_PROMOTER_TYPE_1 = "RNApol_III_promoter_type_1"
PSEUDOGENIC_GENE_SEGMENT = "pseudogenic_gene_segment"
WILD_TYPE_RESCUE_GENE = "wild_type_rescue_gene"
GROUP_1_INTRON_HOMING_ENDONUCLEASE_TARGET_REGION = "group_1_intron_homing_endonuclease_target_region"
RNASE_P_RNA_GENE = "RNase_P_RNA_gene"
V_DJ_CLUSTER = "V_DJ_cluster"
RNASE_P_RNA = "RNase_P_RNA"
CLEAVED_FOR_GPI_ANCHOR_REGION = "cleaved_for_gpi_anchor_region"
ORF = "ORF"
RESCUE_MINI_GENE = "rescue_mini_gene"
HERITABLE_PHENOTYPIC_MARKER = "heritable_phenotypic_marker"
BAC_CLONED_GENOMIC_INSERT = "BAC_cloned_genomic_insert"
POLYPEPTIDE_NICKEL_ION_CONTACT_SITE = "polypeptide_nickel_ion_contact_site"
HAPLOTYPE_BLOCK = "haplotype_block"
THREE_PRIME_INTRON = "three_prime_intron"
INTERMEDIATE_ELEMENT = "intermediate_element"
DESTRUCTION_BOX = "destruction_box"
H3K27_DIMETHYLATION_SITE = "H3K27_dimethylation_site"
S_GNA_OLIGO = "S_GNA_oligo"
GENE_COMPONENT_REGION = "gene_component_region"
V_VJ_CLUSTER = "V_VJ_cluster"
ORIV = "oriV"
SNRNA_GENE = "snRNA_gene"
FOUR_DEMETHYLWYOSINE = "four_demethylwyosine"
TELOMERE = "telomere"
EXTRAMEMBRANE_POLYPEPTIDE_REGION = "extramembrane_polypeptide_region"
RESTRICTION_ENZYME_CLEAVAGE_JUNCTION = "restriction_enzyme_cleavage_junction"
BIDIRECTIONAL_PROMOTER = "bidirectional_promoter"
UPD = "UPD"
SYNTENIC_REGION = "syntenic_region"
ANTIPARALLEL_BETA_STRAND = "antiparallel_beta_strand"
ENHANCER_BOUND_BY_FACTOR = "enhancer_bound_by_factor"
DPE1_MOTIF = "DPE1_motif"
FLOXED_GENE = "floxed_gene"
ALANINE_TRNA_PRIMARY_TRANSCRIPT = "alanine_tRNA_primary_transcript"
NESTED_TRANSPOSON = "nested_transposon"
SHRNA = "shRNA"
RIBOTHYMIDINE = "ribothymidine"
CLONED_GENOMIC_INSERT = "cloned_genomic_insert"
EXPERIMENTAL_FEATURE = "experimental_feature"
PROTEIN_BINDING_SITE = "protein_binding_site"
SCHELLMANN_LOOP_SIX = "schellmann_loop_six"
ENDOSOMAL_LOCALIZATION_SIGNAL = "endosomal_localization_signal"
_8_OXOGUANINE = "8_oxoguanine"
GLYCYL_TRNA = "glycyl_tRNA"
INTROGRESSED_CHROMOSOME_REGION = "introgressed_chromosome_region"
ONE_METHYLADENOSINE = "one_methyladenosine"
FIVE_METHYLAMINOMETHYLURIDINE = "five_methylaminomethyluridine"
BETA_STRAND = "beta_strand"
HOMING_ENDONUCLEASE_BINDING_SITE = "homing_endonuclease_binding_site"
TRYPTOPHANYL_TRNA = "tryptophanyl_tRNA"
CYANELLE_GENE = "cyanelle_gene"
MODIFIED_GUANINE = "modified_guanine"
GALACTOSYL_QUEUOSINE = "galactosyl_queuosine"
MATCH = "match"
POINT_MUTATION = "point_mutation"
HISTONE_4_ACYLATION_SITE = "histone_4_acylation_site"
CDS_PREDICTED = "CDS_predicted"
TRANSCRIPT_BOUND_BY_NUCLEIC_ACID = "transcript_bound_by_nucleic_acid"
BETA_TURN_TYPE_SIX_A_ONE = "beta_turn_type_six_a_one"
TRANSGENIC_TRANSPOSABLE_ELEMENT = "transgenic_transposable_element"
PYROSEQUENCED_READ = "pyrosequenced_read"
ONE_METHYL_THREE_THREE_AMINO_THREE_CARBOXYPROPYL_PSEUDOURIDINE = "one_methyl_three_three_amino_three_carboxypropyl_pseudouridine"
NUCLEOTIDE_MOTIF = "nucleotide_motif"
V_VDJ_J_C_CLUSTER = "V_VDJ_J_C_cluster"
ALANINE = "alanine"
H4K5_ACYLATION_SITE = "H4K5_acylation_site"
RECURSIVE_SPLICE_SITE = "recursive_splice_site"
FOREIGN_GENE = "foreign_gene"
VJ_J_C_CLUSTER = "VJ_J_C_cluster"
FIVE_METHOXYCARBONYLMETHYL_TWO_THIOURIDINE = "five_methoxycarbonylmethyl_two_thiouridine"
H3K9_TRIMETHYLATION_SITE = "H3K9_trimethylation_site"
POLYPEPTIDE_STRUCTURAL_MOTIF = "polypeptide_structural_motif"
MATING_TYPE_REGION_REPLICATION_FORK_BARRIER = "mating_type_region_replication_fork_barrier"
FIVE_FORMYL_TWO_PRIME_O_METHYLCYTIDINE = "five_formyl_two_prime_O_methylcytidine"
COPY_NUMBER_GAIN = "copy_number_gain"
BETA_TURN_TYPE_SIX = "beta_turn_type_six"
RRNA_5_8S = "rRNA_5_8S"
TWO_PRIME_O_METHYLGUANOSINE = "two_prime_O_methylguanosine"
DICISTRONIC_TRANSCRIPT = "dicistronic_transcript"
THREE_PRIME_CODING_EXON_NONCODING_REGION = "three_prime_coding_exon_noncoding_region"
GROUP_IIA_INTRON = "group_IIA_intron"
CAPPED_PRIMARY_TRANSCRIPT = "capped_primary_transcript"
FIVE_PRIME_CODING_EXON_NONCODING_REGION = "five_prime_coding_exon_noncoding_region"
RRNA_CLEAVAGE_SNORNA_PRIMARY_TRANSCRIPT = "rRNA_cleavage_snoRNA_primary_transcript"
A_BOX_TYPE_2 = "A_box_type_2"
PSE_MOTIF = "PSE_motif"
ENGINEERED_GENE = "engineered_gene"
REPEAT_REGION = "repeat_region"
BAC_END = "BAC_end"
LINE_ELEMENT = "LINE_element"
FOREIGN_TRANSPOSABLE_ELEMENT = "foreign_transposable_element"
OCTAMER_MOTIF = "octamer_motif"
PURINE_TO_PYRIMIDINE_TRANSVERSION = "purine_to_pyrimidine_transversion"
VJ_J_CLUSTER = "VJ_J_cluster"
CENTROMERIC_REPEAT = "centromeric_repeat"
MEMBRANE_PEPTIDE_LOOP = "membrane_peptide_loop"
SRP_RNA_PRIMARY_TRANSCRIPT = "SRP_RNA_primary_transcript"
H3K36_TRIMETHYLATION_SITE = "H3K36_trimethylation_site"
CHROMOSOME_PART = "chromosome_part"
GCVB_RNA = "GcvB_RNA"
BETA_TURN_TYPE_SIX_A = "beta_turn_type_six_a"
THREE_METHYLPSEUDOURIDINE = "three_methylpseudouridine"
TILING_PATH = "tiling_path"
PRIMARY_TRANSCRIPT = "primary_transcript"
PROCESSED_PSEUDOGENE = "processed_pseudogene"
FIVE_PRIME_CODING_EXON = "five_prime_coding_exon"
GOLDEN_PATH_FRAGMENT = "golden_path_fragment"
EXPERIMENTAL_RESULT_REGION = "experimental_result_region"
CATMAT_LEFT_HANDED_THREE = "catmat_left_handed_three"
GENOMIC_ISLAND = "genomic_island"
TRANSCRIPTIONAL_CIS_REGULATORY_REGION = "transcriptional_cis_regulatory_region"
QUEUOSINE = "queuosine"
BASE_PAIR = "base_pair"
PRE_MIRNA = "pre_miRNA"
SERINE_THREONINE_MOTIF = "serine_threonine_motif"
RRNA_GENE = "rRNA_gene"
SHINE_DALGARNO_SEQUENCE = "Shine_Dalgarno_sequence"
SILENCED_GENE = "silenced_gene"
VDJ_GENE_SEGMENT = "VDJ_gene_segment"
INVERSION_SITE_PART = "inversion_site_part"
RNA_POLYMERASE_II_TATA_BOX = "RNA_polymerase_II_TATA_box"
SGRNA = "sgRNA"
INTERIOR_INTRON = "interior_intron"
SNV = "SNV"
PREDICTED_GENE = "predicted_gene"
N6_ISOPENTENYLADENOSINE = "N6_isopentenyladenosine"
TERRA = "TERRA"
INTRONIC_REGULATORY_REGION = "intronic_regulatory_region"
ASYMMETRIC_RNA_INTERNAL_LOOP = "asymmetric_RNA_internal_loop"
SL11_ACCEPTOR_SITE = "SL11_acceptor_site"
DSRA_RNA = "DsrA_RNA"
DIF_SITE = "dif_site"
RRNA_18S = "rRNA_18S"
MIR_ENCODING_SNORNA_PRIMARY_TRANSCRIPT = "miR_encoding_snoRNA_primary_transcript"
THREE_PRIME_REPEAT_RECODING_SIGNAL = "three_prime_repeat_recoding_signal"
EDITED_MRNA = "edited_mRNA"
MOBILE_INTRON = "mobile_intron"
TERMINAL_INVERTED_REPEAT_ELEMENT = "terminal_inverted_repeat_element"
ISOLEUCINE = "isoleucine"
TYROSINE = "tyrosine"
MINI_GENE = "mini_gene"
DCAPS_PRIMER = "dCAPS_primer"
DS_RNA_VIRAL_SEQUENCE = "ds_RNA_viral_sequence"
H3K23_DIMETHYLATION_SITE = "H3K23_dimethylation_site"
NO_OUTPUT = "no_output"
TRANSCRIBED_CLUSTER = "transcribed_cluster"
TWO_METHYLADENOSINE = "two_methyladenosine"
N2_2_PRIME_O_DIMETHYLGUANOSINE = "N2_2_prime_O_dimethylguanosine"
RESTRICTION_ENZYME_ASSEMBLY_SCAR = "restriction_enzyme_assembly_scar"
RETINOIC_ACID_RESPONSIVE_ELEMENT = "retinoic_acid_responsive_element"
GENE_SILENCED_BY_DNA_METHYLATION = "gene_silenced_by_DNA_methylation"
OUTRON = "outron"
SNP = "SNP"
TWO_THIO_TWO_PRIME_O_METHYLURIDINE = "two_thio_two_prime_O_methyluridine"
SEQUENCE_REARRANGEMENT_FEATURE = "sequence_rearrangement_feature"
DIHYDROURIDINE = "dihydrouridine"
KOZAK_SEQUENCE = "kozak_sequence"
MOLECULAR_CONTACT_REGION = "molecular_contact_region"
H3K79_TRIMETHYLATION_SITE = "H3K79_trimethylation_site"
FIVE_METHOXYCARBONYLMETHYL_TWO_PRIME_O_METHYLURIDINE = "five_methoxycarbonylmethyl_two_prime_O_methyluridine"
GENE_CASSETTE_ARRAY = "gene_cassette_array"
CAGE_CLUSTER = "CAGE_cluster"
FIVE_PRIME_RECODING_SITE = "five_prime_recoding_site"
VIRAL_SEQUENCE = "viral_sequence"
GLUTAMIC_ACID = "glutamic_acid"
EXON_OF_SINGLE_EXON_GENE = "exon_of_single_exon_gene"
D_LOOP = "D_loop"
MIRNA_PRIMARY_TRANSCRIPT_REGION = "miRNA_primary_transcript_region"
N_TERMINAL_REGION = "n_terminal_region"
MITOCHONDRIAL_SUPERCONTIG = "mitochondrial_supercontig"
N4_2_PRIME_O_DIMETHYLCYTIDINE = "N4_2_prime_O_dimethylcytidine"
CLASS_I_RNA = "class_I_RNA"
CLEAVED_PEPTIDE_REGION = "cleaved_peptide_region"
FIVE_PRIME_RESTRICTION_ENZYME_JUNCTION = "five_prime_restriction_enzyme_junction"
LINEAR_DOUBLE_STRANDED_DNA_CHROMOSOME = "linear_double_stranded_DNA_chromosome"
V_VJ_J_CLUSTER = "V_VJ_J_cluster"
TETRALOOP = "tetraloop"
MOBILE_ELEMENT_INSERTION = "mobile_element_insertion"
REPLICATION_REGULATORY_REGION = "replication_regulatory_region"
D_GENE_RECOMBINATION_FEATURE = "D_gene_recombination_feature"
TWO_PRIME_O_RIBOSYLGUANOSINE_PHOSPHATE = "two_prime_O_ribosylguanosine_phosphate"
C_CLUSTER = "C_cluster"
MICRONUCLEAR_CHROMOSOME = "micronuclear_chromosome"
BACTERIAL_RNAPOL_PROMOTER_SIGMA54 = "bacterial_RNApol_promoter_sigma54"
PEPTIDE_COIL = "peptide_coil"
NC_CONSERVED_REGION = "nc_conserved_region"
RNA_SEQUENCE_SECONDARY_STRUCTURE = "RNA_sequence_secondary_structure"
RNASE_MRP_RNA_GENE = "RNase_MRP_RNA_gene"
REGION = "region"
PLUS_2_TRANSLATIONAL_FRAMESHIFT = "plus_2_translational_frameshift"
J_CLUSTER = "J_cluster"
SEQUENCE_ASSEMBLY = "sequence_assembly"
MRNA = "mRNA"
NONCODING_REGION_OF_EXON = "noncoding_region_of_exon"
DRE_MOTIF = "DRE_motif"
CENTRAL_HYDROPHOBIC_REGION_OF_SIGNAL_PEPTIDE = "central_hydrophobic_region_of_signal_peptide"
BRUNO_RESPONSE_ELEMENT = "Bruno_response_element"
START_CODON = "start_codon"
HSE = "HSE"
D_J_CLUSTER = "D_J_cluster"
TRYPTOPHAN_TRNA_PRIMARY_TRANSCRIPT = "tryptophan_tRNA_primary_transcript"
HOMOLOGOUS_REGION = "homologous_region"
FRT_SITE = "FRT_site"
WOBBLE_BASE_PAIR = "wobble_base_pair"
CLONED_CDNA_INSERT = "cloned_cDNA_insert"
MODIFIED_GUANOSINE = "modified_guanosine"
CLONE_INSERT = "clone_insert"
QTL = "QTL"
ARRET = "ARRET"
PROPEPTIDE = "propeptide"
CODING_CONSERVED_REGION = "coding_conserved_region"
RRNA = "rRNA"
VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_FEATURE = "vertebrate_immune_system_gene_recombination_feature"
CYSTEINE_TRNA_PRIMARY_TRANSCRIPT = "cysteine_tRNA_primary_transcript"
FIVE_TAURINOMETHYLURIDINE = "five_taurinomethyluridine"
G_TO_C_TRANSVERSION = "G_to_C_transversion"
RESCUE_GENE = "rescue_gene"
NDM2_MOTIF = "NDM2_motif"
SONICATE_FRAGMENT = "sonicate_fragment"
ENGINEERED_RESCUE_REGION = "engineered_rescue_region"
INSERTION_BREAKPOINT = "insertion_breakpoint"
GT_DINUCLEOTIDE_REPEAT = "GT_dinucleotide_repeat"
C_BOX = "C_box"
TWO_PRIME_O_METHYLADENOSINE = "two_prime_O_methyladenosine"
H3K56_ACETYLATION_SITE = "H3K56_acetylation_site"
MODIFIED_L_SELENOCYSTEINE = "modified_L_selenocysteine"
FIVE_FORMYLCYTIDINE = "five_formylcytidine"
INSULATOR = "insulator"
BETA_TURN_TYPE_SIX_B = "beta_turn_type_six_b"
POLYPEPTIDE_CATALYTIC_MOTIF = "polypeptide_catalytic_motif"
H3K27_ACYLATION_SITE = "H3K27_acylation_site"
H3R2_DIMETHYLATION_SITE = "H3R2_dimethylation_site"
NUCLEOTIDE_MATCH = "nucleotide_match"
R_LTR_REGION = "R_LTR_region"
SEQUENCE_SECONDARY_STRUCTURE = "sequence_secondary_structure"
INVERSION_SITE = "inversion_site"
SPLICE_REGION = "splice_region"
SCRNA_PRIMARY_TRANSCRIPT = "scRNA_primary_transcript"
TERMINATOR = "terminator"
POLINTON = "polinton"
INTRON_DOMAIN = "intron_domain"
FIVE_PRIME_NONCODING_EXON = "five_prime_noncoding_exon"
TASIRNA = "tasiRNA"
LYSYL_TRNA = "lysyl_tRNA"
LEUCOPLAST_CHROMOSOME = "leucoplast_chromosome"
GENE_WITH_TRANS_SPLICED_TRANSCRIPT = "gene_with_trans_spliced_transcript"
TANDEM_REPEAT = "tandem_repeat"
N2_N2_2_PRIME_O_TRIMETHYLGUANOSINE = "N2_N2_2_prime_O_trimethylguanosine"
ATTC_SITE = "attC_site"
RNAPOL_II_PROMOTER = "RNApol_II_promoter"
MATHEMATICALLY_DEFINED_REPEAT = "mathematically_defined_repeat"
DNA_SEQUENCE_SECONDARY_STRUCTURE = "DNA_sequence_secondary_structure"
SPLICEOSOMAL_INTRON = "spliceosomal_intron"
A_TO_G_TRANSITION = "A_to_G_transition"
H2BK5_MONOMETHYLATION_SITE = "H2BK5_monomethylation_site"
MIR_ENCODING_Y_RNA_PRIMARY_TRANSCRIPT = "miR_encoding_Y_RNA_primary_transcript"
V_D_J_C_CLUSTER = "V_D_J_C_cluster"
GUIDE_RNA = "guide_RNA"
EDITED_TRANSCRIPT_BY_A_TO_I_SUBSTITUTION = "edited_transcript_by_A_to_I_substitution"
IRON_RESPONSIVE_ELEMENT = "iron_responsive_element"
V_SPACER = "V_spacer"
NUCLEAR_MT_PSEUDOGENE = "nuclear_mt_pseudogene"
VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_REARRANGED_GENE_CLUSTER = "vertebrate_immunoglobulin_T_cell_receptor_rearranged_gene_cluster"
CHLOROPLAST_DNA_READ = "chloroplast_DNA_read"
TMRNA_CODING_PIECE = "tmRNA_coding_piece"
BLOCKED_READING_FRAME = "blocked_reading_frame"
SILENT_MATING_TYPE_CASSETTE_ARRAY = "silent_mating_type_cassette_array"
SCHELLMANN_LOOP_SEVEN = "schellmann_loop_seven"
STOP_CODON = "stop_codon"
POLYPEPTIDE_NEST_RIGHT_LEFT_MOTIF = "polypeptide_nest_right_left_motif"
BRED_MOTIF = "BREd_motif"
THREE_PRIME_CLIP = "three_prime_clip"
ORIT = "oriT"
MRNA_CONTIG = "mRNA_contig"
LARGE_SUBUNIT_RRNA = "large_subunit_rRNA"
NUCLEAR_EXPORT_SIGNAL = "nuclear_export_signal"
INVERTED_REPEAT = "inverted_repeat"
LTR_RETROTRANSPOSON = "LTR_retrotransposon"
TELOMERASE_RNA = "telomerase_RNA"
G_TO_T_TRANSVERSION = "G_to_T_transversion"
RH_MAP = "RH_map"
T_LOOP = "T_loop"
MRNA_WITH_FRAMESHIFT = "mRNA_with_frameshift"
SOLO_LTR = "solo_LTR"
BOUNDARY_ELEMENT = "boundary_element"
NUCLEAR_RIM_LOCALIZATION_SIGNAL = "nuclear_rim_localization_signal"
ENGINEERED_TRANSPOSABLE_ELEMENT = "engineered_transposable_element"
ARGINYL_TRNA = "arginyl_tRNA"
TRANSCRIPT = "transcript"
MATING_TYPE_REGION_MOTIF = "mating_type_region_motif"
MODIFIED_L_CYSTEINE = "modified_L_cysteine"
MONOMERIC_REPEAT = "monomeric_repeat"
INTERNAL_TRANSCRIBED_SPACER_REGION = "internal_transcribed_spacer_region"
METHYLINOSINE = "methylinosine"
FLANKING_REGION = "flanking_region"
PI_HELIX = "pi_helix"
ASPARAGINE_TRNA_PRIMARY_TRANSCRIPT = "asparagine_tRNA_primary_transcript"
BINDING_SITE = "binding_site"
X_ELEMENT_COMBINATORIAL_REPEAT = "X_element_combinatorial_repeat"
ALLELICALLY_EXCLUDED_GENE = "allelically_excluded_gene"
MATERNALLY_IMPRINTED_GENE = "maternally_imprinted_gene"
FIVE_AMINOMETHYL_TWO_THIOURIDINE = "five_aminomethyl_two_thiouridine"
MOBILE_GENETIC_ELEMENT = "mobile_genetic_element"
AACCCT_BOX = "AACCCT_box"
VALYL_TRNA = "valyl_tRNA"
EXEMPLAR_MRNA = "exemplar_mRNA"
ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT_GENE = "engineered_foreign_transposable_element_gene"
TRANSITION = "transition"
DISPERSED_REPEAT = "dispersed_repeat"
CASSETTE_PSEUDOGENE = "cassette_pseudogene"
EXONIC_SPLICE_ENHANCER = "exonic_splice_enhancer"
PYRROLYSINE = "pyrrolysine"
FIVE_METHOXYURIDINE = "five_methoxyuridine"
MRNA_WITH_MINUS_1_FRAMESHIFT = "mRNA_with_minus_1_frameshift"
H3K9_ACETYLATION_SITE = "H3K9_acetylation_site"
ANDROGEN_RESPONSE_ELEMENT = "androgen_response_element"
APTAMER = "aptamer"
FIVE_PRIME_D_NONAMER = "five_prime_D_nonamer"
GENE_SILENCED_BY_HISTONE_METHYLATION = "gene_silenced_by_histone_methylation"
OPERATOR = "operator"
DCE = "DCE"
CONSENSUS_GDNA = "consensus_gDNA"
CDS_SUPPORTED_BY_EST_OR_CDNA_DATA = "CDS_supported_by_EST_or_cDNA_data"
FIVE_TAURINOMETHYL_TWO_THIOURIDINE = "five_taurinomethyl_two_thiouridine"
SNORNA_PRIMARY_TRANSCRIPT = "snoRNA_primary_transcript"
RDNA_REPLICATION_FORK_BARRIER = "rDNA_replication_fork_barrier"
VDJ_J_C_CLUSTER = "VDJ_J_C_cluster"
H3R2_MONOMETHYLATION_SITE = "H3R2_monomethylation_site"
STOP_CODON_REDEFINED_AS_PYRROLYSINE = "stop_codon_redefined_as_pyrrolysine"
EDITED_TRANSCRIPT_FEATURE = "edited_transcript_feature"
SL4_ACCEPTOR_SITE = "SL4_acceptor_site"
THREE_THREE_AMINO_THREE_CARBOXYPROPYL_URIDINE = "three_three_amino_three_carboxypropyl_uridine"
CATMAT_RIGHT_HANDED_FOUR = "catmat_right_handed_four"
U5_FIVE_PRIME_LTR_REGION = "U5_five_prime_LTR_region"
TRNA_PRIMARY_TRANSCRIPT = "tRNA_primary_transcript"
DNA_BINDING_SITE = "DNA_binding_site"
FINGERPRINT_MAP = "fingerprint_map"
GENE_WITH_DICISTRONIC_TRANSCRIPT = "gene_with_dicistronic_transcript"
UGA_STOP_CODON_SIGNAL = "UGA_stop_codon_signal"
H3K23_ACYLATION_SITE = "H3K23_acylation_site"
THREE_PRIME_RECODING_SITE = "three_prime_recoding_site"
INTERIOR_EXON = "interior_exon"
THREE_PRIME_UTR_INTRON = "three_prime_UTR_intron"
MAXICIRCLE_GENE = "maxicircle_gene"
N2_N2_DIMETHYLGUANOSINE = "N2_N2_dimethylguanosine"
T3_RNA_POLYMERASE_PROMOTER = "T3_RNA_Polymerase_Promoter"
J_SPACER = "J_spacer"
SITE_SPECIFIC_RECOMBINATION_TARGET_REGION = "site_specific_recombination_target_region"
RECODING_PSEUDOKNOT = "recoding_pseudoknot"
GENOMIC_DNA_READ = "genomic_DNA_read"
MICF_RNA = "MicF_RNA"
IMMUNOGLOBULIN_REGION = "immunoglobulin_region"
SPLICEOSOMAL_INTRON_REGION = "spliceosomal_intron_region"
MRNA_RECODED_BY_CODON_REDEFINITION = "mRNA_recoded_by_codon_redefinition"
RASIRNA = "rasiRNA"
TANDEM_DUPLICATION = "tandem_duplication"
POLYPEPTIDE_COBALT_ION_CONTACT_SITE = "polypeptide_cobalt_ion_contact_site"
INTRINSICALLY_UNSTRUCTURED_POLYPEPTIDE_REGION = "intrinsically_unstructured_polypeptide_region"
CENTROMERE_DNA_ELEMENT_II = "centromere_DNA_Element_II"
SEQUENCE_ALTERATION = "sequence_alteration"
VECTOR_REPLICON = "vector_replicon"
FIVE_PRIME_EST = "five_prime_EST"
AMPLIFICATION_ORIGIN = "amplification_origin"
CONSENSUS_REGION = "consensus_region"
N6_N6_2_PRIME_O_TRIMETHYLADENOSINE = "N6_N6_2_prime_O_trimethyladenosine"
A_TO_C_TRANSVERSION = "A_to_C_transversion"
CHROMOPLAST_GENE = "chromoplast_gene"
COILED_COIL = "coiled_coil"
TRYPTOPHAN = "tryptophan"
V_D_DJ_C_CLUSTER = "V_D_DJ_C_cluster"
ASX_TURN_LEFT_HANDED_TYPE_TWO = "asx_turn_left_handed_type_two"
X_REGION = "X_region"
TRANSPOSON_FRAGMENT = "transposon_fragment"
NONCODING_EXON = "noncoding_exon"
SL2_ACCEPTOR_SITE = "SL2_acceptor_site"
HISTONE_MODIFICATION = "histone_modification"
U_BOX = "U_box"
DOUBLE_STRANDED_RNA_CHROMOSOME = "double_stranded_RNA_chromosome"
PYRROLYSYL_TRNA = "pyrrolysyl_tRNA"
NUCLEASE_BINDING_SITE = "nuclease_binding_site"
MIRNA = "miRNA"
FIVE_PRIME_FLANKING_REGION = "five_prime_flanking_region"
CRISPR = "CRISPR"
SARCIN_LIKE_RNA_MOTIF = "sarcin_like_RNA_motif"
RHO_DEPENDENT_BACTERIAL_TERMINATOR = "rho_dependent_bacterial_terminator"
H3K4_MONOMETHYLATION_SITE = "H3K4_monomethylation_site"
RNA_POLYMERASE_PROMOTER = "RNA_polymerase_promoter"
SINGLE_STRANDED_RNA_CHROMOSOME = "single_stranded_RNA_chromosome"
N4_METHYLCYTIDINE = "N4_methylcytidine"
THREE_PRIME_RESTRICTION_ENZYME_JUNCTION = "three_prime_restriction_enzyme_junction"
POLYPEPTIDE_TURN_MOTIF = "polypeptide_turn_motif"
THREE_PRIME_UST = "three_prime_UST"
DCE_SI = "DCE_SI"
ST_TURN_LEFT_HANDED_TYPE_ONE = "st_turn_left_handed_type_one"
HIGH_IDENTITY_REGION = "high_identity_region"
H3K20_TRIMETHYLATION_SITE = "H3K20_trimethylation_site"
MODIFIED_ADENINE = "modified_adenine"
SEQUENCE_CONFLICT = "sequence_conflict"
AMINO_ACID = "amino_acid"
ANTICODON = "anticodon"
SIMPLE_SEQUENCE_LENGTH_VARIATION = "simple_sequence_length_variation"
BAC_READ_CONTIG = "BAC_read_contig"
SIGNAL_ANCHOR = "signal_anchor"
INTERNAL_RIBOSOME_ENTRY_SITE = "internal_ribosome_entry_site"
V_D_DJ_J_CLUSTER = "V_D_DJ_J_cluster"
SRP_RNA_GENE = "SRP_RNA_gene"
VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_SEGMENT = "vertebrate_immunoglobulin_T_cell_receptor_segment"
INDUCIBLE_PROMOTER = "inducible_promoter"
V_DJ_J_CLUSTER = "V_DJ_J_cluster"
DUPLICATED_PSEUDOGENE = "duplicated_pseudogene"
_5_HYDROXYMETHYLCYTOSINE = "5_hydroxymethylcytosine"
REGIONAL_CENTROMERE_INNER_REPEAT_REGION = "regional_centromere_inner_repeat_region"
FIVE_PRIME_LTR = "five_prime_LTR"
U4ATAC_SNRNA = "U4atac_snRNA"
MATING_TYPE_M_BOX = "mating_type_M_box"
EXONIC_SPLICING_SILENCER = "exonic_splicing_silencer"
ANCHOR_REGION = "anchor_region"
CONSERVED_REGION = "conserved_region"
MAXICIRCLE = "maxicircle"
TNAORF = "tnaORF"
GENE_WITH_POLYADENYLATED_MRNA = "gene_with_polyadenylated_mRNA"
SMALL_SUBUNIT_RRNA = "small_subunit_rRNA"
FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = "five_prime_sticky_end_restriction_enzyme_cleavage_site"
SL8_ACCEPTOR_SITE = "SL8_acceptor_site"
REMARK = "remark"
SCRNA_GENE = "scRNA_gene"
V_CLUSTER = "V_cluster"
INTRON = "intron"
INTEIN_ENCODING_REGION = "intein_encoding_region"
OVERLAPPING_EST_SET = "overlapping_EST_set"
CURE = "CuRE"
GENE_SILENCED_BY_HISTONE_MODIFICATION = "gene_silenced_by_histone_modification"
URIDINE_FIVE_OXYACETIC_ACID_METHYL_ESTER = "uridine_five_oxyacetic_acid_methyl_ester"
LYSINE_TRNA_PRIMARY_TRANSCRIPT = "lysine_tRNA_primary_transcript"
TRANS_SPLICED_MRNA = "trans_spliced_mRNA"
RHO_INDEPENDENT_BACTERIAL_TERMINATOR = "rho_independent_bacterial_terminator"
MODIFIED_L_ARGININE = "modified_L_arginine"
GENE_ARRAY = "gene_array"
MODIFIED_L_ISOLEUCINE = "modified_L_isoleucine"
T_TO_C_TRANSITION = "T_to_C_transition"
N4_ACETYL_2_PRIME_O_METHYLCYTIDINE = "N4_acetyl_2_prime_O_methylcytidine"
TWO_PRIME_O_METHYLCYTIDINE = "two_prime_O_methylcytidine"
TRNA = "tRNA"
VALINE = "valine"
DELETION_JUNCTION = "deletion_junction"
PSEUDOGENIC_TRANSCRIPT = "pseudogenic_transcript"
COPY_NUMBER_LOSS = "copy_number_loss"
U4_SNRNA = "U4_snRNA"
VALINE_TRNA_PRIMARY_TRANSCRIPT = "valine_tRNA_primary_transcript"
COSMID = "cosmid"
THREE_TEN_HELIX = "three_ten_helix"
EUKARYOTIC_TERMINATOR = "eukaryotic_terminator"
FIVE_METHYL_2_THIOURIDINE = "five_methyl_2_thiouridine"
FIVE_CARBOXYHYDROXYMETHYL_URIDINE_METHYL_ESTER = "five_carboxyhydroxymethyl_uridine_methyl_ester"
CONJUGATIVE_TRANSPOSON = "conjugative_transposon"
PRIMER = "primer"
EPISOME = "episome"
ULTRACONTIG = "ultracontig"
DJ_J_C_CLUSTER = "DJ_J_C_cluster"
COMPLEX_STRUCTURAL_ALTERATION = "complex_structural_alteration"
POLYPEPTIDE_ZINC_ION_CONTACT_SITE = "polypeptide_zinc_ion_contact_site"
CODING_END = "coding_end"
BRANCH_SITE = "branch_site"
TWO_PRIME_O_METHYLINOSINE = "two_prime_O_methylinosine"
METHIONINE = "methionine"
UTR_INTRON = "UTR_intron"
PROMOTER_ELEMENT = "promoter_element"
POST_TRANSLATIONALLY_MODIFIED_REGION = "post_translationally_modified_region"
D_DJ_J_CLUSTER = "D_DJ_J_cluster"
FOLDBACK_ELEMENT = "foldback_element"
G_TO_A_TRANSITION = "G_to_A_transition"
GLUTAMINYL_TRNA = "glutaminyl_tRNA"
FIVE_CARBOXYMETHYLURIDINE = "five_carboxymethyluridine"
GENE_WITH_DICISTRONIC_MRNA = "gene_with_dicistronic_mRNA"
ATTP_SITE = "attP_site"
TRANSLOCATION_BREAKPOINT = "translocation_breakpoint"
REGULATORY_PROMOTER_ELEMENT = "regulatory_promoter_element"
MATERNAL_UNIPARENTAL_DISOMY = "maternal_uniparental_disomy"
K_TURN_RNA_MOTIF = "K_turn_RNA_motif"
NON_CYTOPLASMIC_POLYPEPTIDE_REGION = "non_cytoplasmic_polypeptide_region"
CIRCULAR_DOUBLE_STRANDED_RNA_CHROMOSOME = "circular_double_stranded_RNA_chromosome"
GAMMA_TURN_CLASSIC = "gamma_turn_classic"
RNAPOL_III_PROMOTER = "RNApol_III_promoter"
PSEUDOGENIC_EXON = "pseudogenic_exon"
NATURAL_PLASMID = "natural_plasmid"
APICOPLAST_GENE = "apicoplast_gene"
EDITING_DOMAIN = "editing_domain"
INTRON_BASE_5 = "intron_base_5"
MODIFIED_RNA_BASE_FEATURE = "modified_RNA_base_feature"
LEUCINE = "leucine"
TRANSGENE = "transgene"
THREE_PRIME_FLANKING_REGION = "three_prime_flanking_region"
H3K4_METHYLATION_SITE = "H3K4_methylation_site"
INTERCHROMOSOMAL_BREAKPOINT = "interchromosomal_breakpoint"
PHOSPHORYLATION_SITE = "phosphorylation_site"
REGIONAL_CENTROMERE = "regional_centromere"
RAPD = "RAPD"
H3K14_ACETYLATION_SITE = "H3K14_acetylation_site"
POLYPEPTIDE_DOMAIN = "polypeptide_domain"
HOOGSTEEN_BASE_PAIR = "Hoogsteen_base_pair"
PROXIMAL_PROMOTER_ELEMENT = "proximal_promoter_element"
C_D_BOX_SNORNA = "C_D_box_snoRNA"
RESTRICTION_ENZYME_RECOGNITION_SITE = "restriction_enzyme_recognition_site"
POLYPEPTIDE_MOTIF = "polypeptide_motif"
CSRB_RSMB_RNA = "CsrB_RsmB_RNA"
TMRNA = "tmRNA"
POSITIONAL_CANDIDATE_GENE = "positional_candidate_gene"
LARIAT_INTRON = "lariat_intron"
LYSINE = "lysine"
ATTR_SITE = "attR_site"
SUGAR_EDGE_BASE_PAIR = "sugar_edge_base_pair"
BETA_TURN_LEFT_HANDED_TYPE_ONE = "beta_turn_left_handed_type_one"
CHIP_SEQ_REGION = "ChIP_seq_region"
TRANSGENIC_INSERTION = "transgenic_insertion"
SIGNAL_PEPTIDE = "signal_peptide"
W_REGION = "W_region"
S_REGION = "S_region"
RRNA_21S = "rRNA_21S"
FUSION_GENE = "fusion_gene"
CHROMOSOMAL_STRUCTURAL_ELEMENT = "chromosomal_structural_element"
N6_HYDROXYNORVALYLCARBAMOYLADENOSINE = "N6_hydroxynorvalylcarbamoyladenosine"
H2B_UBIQUITINATION_SITE = "H2B_ubiquitination_site"
RST_MATCH = "RST_match"
CCA_TAIL = "CCA_tail"
SCHELLMANN_LOOP = "schellmann_loop"
GENE_SILENCED_BY_DNA_MODIFICATION = "gene_silenced_by_DNA_modification"
ENGINEERED_INSERT = "engineered_insert"
UNITARY_PSEUDOGENE = "unitary_pseudogene"
FIVE_PRIME_RST = "five_prime_RST"
U12_SNRNA = "U12_snRNA"
U14_SNORNA = "U14_snoRNA"
MODIFIED_L_LEUCINE = "modified_L_leucine"
RESOLUTION_SITE = "resolution_site"
SL7_ACCEPTOR_SITE = "SL7_acceptor_site"
PROBE = "probe"
POLYPEPTIDE_SECONDARY_STRUCTURE = "polypeptide_secondary_structure"
POLYPEPTIDE_SEQUENCING_INFORMATION = "polypeptide_sequencing_information"
V_GENE_RECOMBINATION_FEATURE = "V_gene_recombination_feature"
FIVE_CARBAMOYLMETHYL_TWO_PRIME_O_METHYLURIDINE = "five_carbamoylmethyl_two_prime_O_methyluridine"
TRINUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE = "trinucleotide_repeat_microsatellite_feature"
CLONED_REGION = "cloned_region"
MIR_ENCODING_SHRNA_PRIMARY_TRANSCRIPT = "miR_encoding_shRNA_primary_transcript"
GENE_SILENCED_BY_HISTONE_DEACETYLATION = "gene_silenced_by_histone_deacetylation"
CIRCULAR_SINGLE_STRANDED_RNA_CHROMOSOME = "circular_single_stranded_RNA_chromosome"
H3K9_MONOMETHYLATION_SITE = "H3K9_monomethylation_site"
ENDOGENOUS_RETROVIRAL_GENE = "endogenous_retroviral_gene"
FIVE_METHYLCYTIDINE = "five_methylcytidine"
STEM_LOOP = "stem_loop"
PROLYL_TRNA = "prolyl_tRNA"
SCRNA = "scRNA"
ORTHOLOGOUS_REGION = "orthologous_region"
FOSMID = "fosmid"
INTEIN = "intein"
GENE_WITH_MRNA_RECODED_BY_TRANSLATIONAL_BYPASS = "gene_with_mRNA_recoded_by_translational_bypass"
IRON_REPRESSED_GATA_ELEMENT = "iron_repressed_GATA_element"
HISTIDINE = "histidine"
POLYPEPTIDE_MAGNESIUM_ION_CONTACT_SITE = "polypeptide_magnesium_ion_contact_site"
OKAZAKI_FRAGMENT = "Okazaki_fragment"
UST = "UST"
H4K12_ACETYLATION_SITE = "H4K12_acetylation_site"
POLYPEPTIDE_CONSERVED_MOTIF = "polypeptide_conserved_motif"
NATURAL_VARIANT_SITE = "natural_variant_site"
H3K4_DIMETHYLATION_SITE = "H3K4_dimethylation_site"
OVERLAPPING_FEATURE_SET = "overlapping_feature_set"
GLYCINE_TRNA_PRIMARY_TRANSCRIPT = "glycine_tRNA_primary_transcript"
RIBONUCLEASE_SITE = "ribonuclease_site"
PLASTID_GENE = "plastid_gene"
VD_GENE_SEGMENT = "VD_gene_segment"
V_VDJ_J_CLUSTER = "V_VDJ_J_cluster"
SINGLE_STRANDED_DNA_CHROMOSOME = "single_stranded_DNA_chromosome"
POLYPEPTIDE_NEST_LEFT_RIGHT_MOTIF = "polypeptide_nest_left_right_motif"
ASX_MOTIF = "asx_motif"
RNAPOL_I_PROMOTER = "RNApol_I_promoter"
ASSEMBLY_ERROR_CORRECTION = "assembly_error_correction"
PHENYLALANYL_TRNA = "phenylalanyl_tRNA"
BLUNT_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = "blunt_end_restriction_enzyme_cleavage_site"
UAG_STOP_CODON_SIGNAL = "UAG_stop_codon_signal"
POLYPEPTIDE_REGION = "polypeptide_region"
RESTRICTION_ENZYME_THREE_PRIME_SINGLE_STRAND_OVERHANG = "restriction_enzyme_three_prime_single_strand_overhang"
ST_TURN_RIGHT_HANDED_TYPE_TWO = "st_turn_right_handed_type_two"
RNA_INTERNAL_LOOP = "RNA_internal_loop"
POLYA_SIGNAL_SEQUENCE = "polyA_signal_sequence"
PCR_PRODUCT = "PCR_product"
IDNA = "iDNA"
D_CLUSTER = "D_cluster"
POLYA_SITE = "polyA_site"
ISRE = "ISRE"
MONOCISTRONIC_TRANSCRIPT = "monocistronic_transcript"
DHU_LOOP = "DHU_loop"
THREE_PRIME_D_HEPTAMER = "three_prime_D_heptamer"
TASIRNA_PRIMARY_TRANSCRIPT = "tasiRNA_primary_transcript"
ZINC_FINGER_BINDING_SITE = "zinc_finger_binding_site"
REPEAT_COMPONENT = "repeat_component"
LINEAR_DOUBLE_STRANDED_RNA_CHROMOSOME = "linear_double_stranded_RNA_chromosome"
URIDINE_FIVE_OXYACETIC_ACID = "uridine_five_oxyacetic_acid"
C_TO_G_TRANSVERSION = "C_to_G_transversion"
ASPARTIC_ACID_TRNA_PRIMARY_TRANSCRIPT = "aspartic_acid_tRNA_primary_transcript"
CTCF_BINDING_SITE = "CTCF_binding_site"
ANTICODON_LOOP = "anticodon_loop"
FOUR_THIOURIDINE = "four_thiouridine"
POLYMORPHIC_PSEUDOGENE = "polymorphic_pseudogene"
R_FIVE_PRIME_LTR_REGION = "R_five_prime_LTR_region"
POLYA_PRIMED_CDNA_CLONE = "polyA_primed_cDNA_clone"
SPLICE_SITE = "splice_site"
DMV2_MOTIF = "DMv2_motif"
Z2_REGION = "Z2_region"
FIVE_PRIME_CODING_EXON_CODING_REGION = "five_prime_coding_exon_coding_region"
FIVE_CARBAMOYLMETHYLURIDINE = "five_carbamoylmethyluridine"
FIVE_PRIME_CLIP = "five_prime_clip"
CRYPTIC_GENE = "cryptic_gene"
SEVEN_METHYLGUANINE = "seven_methylguanine"
RETROGENE = "retrogene"
H3K27_METHYLATION_SITE = "H3K27_methylation_site"
NC_PRIMARY_TRANSCRIPT = "nc_primary_transcript"
GROUP_IIB_INTRON = "group_IIB_intron"
MODIFIED_L_VALINE = "modified_L_valine"
H3K79_METHYLATION_SITE = "H3K79_methylation_site"
PLASMID_GENE = "plasmid_gene"
GENE_SILENCED_BY_RNA_INTERFERENCE = "gene_silenced_by_RNA_interference"
CIS_REGULATORY_FRAMESHIFT_ELEMENT = "cis_regulatory_frameshift_element"
REGULATORY_REGION = "regulatory_region"
THREE_PRIME_D_RECOMBINATION_SIGNAL_SEQUENCE = "three_prime_D_recombination_signal_sequence"
MINUS_12_SIGNAL = "minus_12_signal"
POSITIVE_SENSE_SSRNA_VIRAL_SEQUENCE = "positive_sense_ssRNA_viral_sequence"
RECOMBINATION_ENHANCER = "recombination_enhancer"
RNASE_MRP_RNA = "RNase_MRP_RNA"
GENE_WITH_MRNA_WITH_FRAMESHIFT = "gene_with_mRNA_with_frameshift"
ASX_TURN_RIGHT_HANDED_TYPE_TWO = "asx_turn_right_handed_type_two"
ENHANCERRNA = "enhancerRNA"
CHLOROPLAST_CHROMOSOME = "chloroplast_chromosome"
FORKHEAD_MOTIF = "forkhead_motif"
ATTCTN_SITE = "attCtn_site"
PHAGE_RNA_POLYMERASE_PROMOTER = "Phage_RNA_Polymerase_Promoter"
PSEUDOKNOT = "pseudoknot"
FIVE_METHYLDIHYDROURIDINE = "five_methyldihydrouridine"
ATTI_SITE = "attI_site"
MODIFIED_CYTIDINE = "modified_cytidine"
INDEL = "indel"
D_GENE_SEGMENT = "D_gene_segment"
CRYPTIC_SPLICE_SITE = "cryptic_splice_site"
ENGINEERED_REGION = "engineered_region"
SYMBIOSIS_ISLAND = "symbiosis_island"
MIRNA_GENE = "miRNA_gene"
ALANYL_TRNA = "alanyl_tRNA"
DYE_TERMINATOR_READ = "dye_terminator_read"
FIVE_TWO_PRIME_O_DIMETHYLURIDINE = "five_two_prime_O_dimethyluridine"
FIVE_METHYLAMINOMETHYL_TWO_THIOURIDINE = "five_methylaminomethyl_two_thiouridine"
TAG = "tag"
COMPOSITIONALLY_BIASED_REGION_OF_PEPTIDE = "compositionally_biased_region_of_peptide"
RRNA_CLEAVAGE_RNA = "rRNA_cleavage_RNA"
POLYPEPTIDE_IRON_ION_CONTACT_SITE = "polypeptide_iron_ion_contact_site"
INTERNAL_GUIDE_SEQUENCE = "internal_guide_sequence"
J_C_CLUSTER = "J_C_cluster"
ENGINEERED_FOREIGN_REPETITIVE_ELEMENT = "engineered_foreign_repetitive_element"
ASPARTYL_TRNA = "aspartyl_tRNA"
FRAGMENT_ASSEMBLY = "fragment_assembly"
CONTIG_COLLECTION = "contig_collection"
MNP = "MNP"
GROUP_III_INTRON = "group_III_intron"
CLONE_END = "clone_end"
RNAPOL_III_PROMOTER_TYPE_3 = "RNApol_III_promoter_type_3"
H3K36_ACETYLATION_SITE = "H3K36_acetylation_site"
PSEUDOURIDYLATION_GUIDE_SNORNA = "pseudouridylation_guide_snoRNA"
PROTEIN_HMM_MATCH = "protein_hmm_match"
LYSIDINE = "lysidine"
CANONICAL_THREE_PRIME_SPLICE_SITE = "canonical_three_prime_splice_site"
ATTL_SITE = "attL_site"
SIRNA = "siRNA"
ALTERNATIVELY_SPLICED_TRANSCRIPT = "alternatively_spliced_transcript"
MODIFIED_L_TRYPTOPHAN = "modified_L_tryptophan"
SL9_ACCEPTOR_SITE = "SL9_acceptor_site"
REPLICON = "replicon"
MODIFIED_L_TYROSINE = "modified_L_tyrosine"
DATABANK_ENTRY = "databank_entry"
LEFT_HANDED_PEPTIDE_HELIX = "left_handed_peptide_helix"
CLONE_INSERT_START = "clone_insert_start"
MEMBRANE_STRUCTURE = "membrane_structure"
SNORNA = "snoRNA"
REARRANGEMENT_BREAKPOINT = "rearrangement_breakpoint"
C_TO_A_TRANSVERSION = "C_to_A_transversion"
H4K91_ACETYLATION_SITE = "H4K91_acetylation_site"
MODIFIED_L_PROLINE = "modified_L_proline"
DECAYED_EXON = "decayed_exon"
EXTENDED_INTRONIC_SPLICE_REGION = "extended_intronic_splice_region"
V_DJ_C_CLUSTER = "V_DJ_C_cluster"
GENE_WITH_DICISTRONIC_PRIMARY_TRANSCRIPT = "gene_with_dicistronic_primary_transcript"
TOPOLOGICALLY_DEFINED_REGION = "topologically_defined_region"
U5_SNRNA = "U5_snRNA"
N2_METHYLGUANOSINE = "N2_methylguanosine"
PROTEIN_MATCH = "protein_match"
PLUS_1_TRANSLATIONAL_FRAMESHIFT = "plus_1_translational_frameshift"
PIP_BOX = "PIP_box"
RRNA_SMALL_SUBUNIT_PRIMARY_TRANSCRIPT = "rRNA_small_subunit_primary_transcript"
METHYLATED_ADENINE = "methylated_adenine"
DUPLICATION = "duplication"
RRNA_5S = "rRNA_5S"
H3K4_TRIMETHYLATION = "H3K4_trimethylation"
INTERIOR_CODING_EXON = "interior_coding_exon"
CDS_REGION = "CDS_region"
RETRON = "retron"
U12_INTRON = "U12_intron"
RIGHT_HANDED_PEPTIDE_HELIX = "right_handed_peptide_helix"
RNAPOL_II_CORE_PROMOTER = "RNApol_II_core_promoter"
SINE_ELEMENT = "SINE_element"
EXTENDED_CIS_SPLICE_SITE = "extended_cis_splice_site"
CRYPTIC_PROPHAGE = "cryptic_prophage"
METHYLATED_CYTOSINE = "methylated_cytosine"
STREP_MOTIF = "STREP_motif"
UNTRANSLATED_REGION_POLYCISTRONIC_MRNA = "untranslated_region_polycistronic_mRNA"
EXTRACHROMOSOMAL_MOBILE_GENETIC_ELEMENT = "extrachromosomal_mobile_genetic_element"
UTR_REGION = "UTR_region"
HYDROXYWYBUTOSINE = "hydroxywybutosine"
MODIFIED_GLYCINE = "modified_glycine"
MITOCHONDRIAL_CONTIG = "mitochondrial_contig"
PARALLEL_BETA_STRAND = "parallel_beta_strand"
RIBOSWITCH = "riboswitch"
NON_PROCESSED_PSEUDOGENE = "non_processed_pseudogene"
ENGINEERED_FOREIGN_REGION = "engineered_foreign_region"
DIRECT_TANDEM_DUPLICATION = "direct_tandem_duplication"
STOP_CODON_REDEFINED_AS_SELENOCYSTEINE = "stop_codon_redefined_as_selenocysteine"
MICROSATELLITE = "microsatellite"
TELOMERIC_TRANSCRIPT = "telomeric_transcript"
T_TO_G_TRANSVERSION = "T_to_G_transversion"
Y_PRIME_ELEMENT = "Y_prime_element"
MIRNA_ANTIGUIDE = "miRNA_antiguide"
MINICIRCLE = "minicircle"
SL12_ACCEPTOR_SITE = "SL12_acceptor_site"
TEMPLATE_REGION = "template_region"
DNA_APTAMER = "DNA_aptamer"
MITOCHONDRIAL_DNA_READ = "mitochondrial_DNA_read"
POSSIBLE_BASE_CALL_ERROR = "possible_base_call_error"
PEPTIDE_HELIX = "peptide_helix"
FIVE_METHOXYCARBONYLMETHYLURIDINE = "five_methoxycarbonylmethyluridine"
THREE_PRIME_TERMINAL_INVERTED_REPEAT = "three_prime_terminal_inverted_repeat"
MODIFIED_L_LYSINE = "modified_L_lysine"
TCS_ELEMENT = "TCS_element"
ANTISENSE_LNCRNA = "antisense_lncRNA"
U3_FIVE_PRIME_LTR_REGION = "U3_five_prime_LTR_region"
PATERNAL_UNIPARENTAL_DISOMY = "paternal_uniparental_disomy"
PURINE_TRANSITION = "purine_transition"
TRANSPOSABLE_ELEMENT_GENE = "transposable_element_gene"
MRNA_REGION = "mRNA_region"
CLEAVED_INITIATOR_METHIONINE = "cleaved_initiator_methionine"
RRNA_16S = "rRNA_16S"
EPOXYQUEUOSINE = "epoxyqueuosine"
MIRNA_PRIMARY_TRANSCRIPT = "miRNA_primary_transcript"
INTEGRON = "integron"
THREE_METHYLCYTIDINE = "three_methylcytidine"
EDITING_BLOCK = "editing_block"
CHROMOSOME_BAND = "chromosome_band"
CODON = "codon"
ASX_TURN_LEFT_HANDED_TYPE_ONE = "asx_turn_left_handed_type_one"
METHYLATION_GUIDE_SNORNA = "methylation_guide_snoRNA"
TRANSCRIPTION_END_SITE = "transcription_end_site"
SL1_ACCEPTOR_SITE = "SL1_acceptor_site"
U1_SNRNA = "U1_snRNA"
CONSTITUTIVE_PROMOTER = "constitutive_promoter"
RECOMBINATIONALLY_INVERTED_GENE = "recombinationally_inverted_gene"
DCE_SIII = "DCE_SIII"
RRNA_25S = "rRNA_25S"
MONOCISTRONIC_PRIMARY_TRANSCRIPT = "monocistronic_primary_transcript"
G_QUARTET = "G_quartet"
LOXP_SITE = "loxP_site"
FIVE_CARBOXYMETHYLAMINOMETHYLURIDINE = "five_carboxymethylaminomethyluridine"
GENOMICALLY_CONTAMINATED_CDNA_CLONE = "genomically_contaminated_cDNA_clone"
ISOLEUCINE_TRNA_PRIMARY_TRANSCRIPT = "isoleucine_tRNA_primary_transcript"
RNAPOL_III_PROMOTER_TYPE_2 = "RNApol_III_promoter_type_2"
MODIFIED_L_ASPARTIC_ACID = "modified_L_aspartic_acid"
PROMOTER_TARGETING_SEQUENCE = "promoter_targeting_sequence"
TWO_METHYLTHIO_N6_ISOPENTENYLADENOSINE = "two_methylthio_N6_isopentenyladenosine"
LINEAR_SINGLE_STRANDED_DNA_CHROMOSOME = "linear_single_stranded_DNA_chromosome"
Y_REGION = "Y_region"
SPLICE_ENHANCER = "splice_enhancer"
SMALL_REGULATORY_NCRNA = "small_regulatory_ncRNA"
PIRNA = "piRNA"
RDNA_INTERGENIC_SPACER_ELEMENT = "rDNA_intergenic_spacer_element"
DNA_TRANSPOSON = "DNA_transposon"
CODING_REGION_OF_EXON = "coding_region_of_exon"
ONE_METHYLINOSINE = "one_methylinosine"
DJ_C_CLUSTER = "DJ_C_cluster"
THREE_PRIME_CODING_EXON = "three_prime_coding_exon"
PYRIMIDINE_TO_PURINE_TRANSVERSION = "pyrimidine_to_purine_transversion"
TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT = "transcript_with_translational_frameshift"
METHYLATION_GUIDE_SNORNA_PRIMARY_TRANSCRIPT = "methylation_guide_snoRNA_primary_transcript"
DICISTRONIC_MRNA = "dicistronic_mRNA"
REAGENT = "reagent"
PHENYLALANINE_TRNA_PRIMARY_TRANSCRIPT = "phenylalanine_tRNA_primary_transcript"
HISTONE_ACETYLATION_SITE = "histone_acetylation_site"
NON_LTR_RETROTRANSPOSON = "non_LTR_retrotransposon"
ALTERNATE_SEQUENCE_SITE = "alternate_sequence_site"
N2_N2_7_TRIMETHYLGUANOSINE = "N2_N2_7_trimethylguanosine"
FIVE_METHYLURIDINE = "five_methyluridine"
TSS = "TSS"
FIVE_METHYLAMINOMETHYL_TWO_SELENOURIDINE = "five_methylaminomethyl_two_selenouridine"
ANTISENSE_PRIMARY_TRANSCRIPT = "antisense_primary_transcript"
TETRANUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE = "tetranucleotide_repeat_microsatellite_feature"
RESTRICTION_ENZYME_FIVE_PRIME_SINGLE_STRAND_OVERHANG = "restriction_enzyme_five_prime_single_strand_overhang"
CANONICAL_FIVE_PRIME_SPLICE_SITE = "canonical_five_prime_splice_site"
ENHANCER = "enhancer"
POLYPEPTIDE_CONSERVED_REGION = "polypeptide_conserved_region"
TRANSLATION_REGULATORY_REGION = "translation_regulatory_region"
ENGINEERED_FUSION_GENE = "engineered_fusion_gene"
GENE_WITH_POLYCISTRONIC_TRANSCRIPT = "gene_with_polycistronic_transcript"
PHAGE_SEQUENCE = "phage_sequence"
MTE = "MTE"
CPG_ISLAND = "CpG_island"
POLYPEPTIDE_VARIATION_SITE = "polypeptide_variation_site"
KEN_BOX = "KEN_box"
FIVE_CARBOXYMETHYLAMINOMETHYL_TWO_THIOURIDINE = "five_carboxymethylaminomethyl_two_thiouridine"
SL5_ACCEPTOR_SITE = "SL5_acceptor_site"
POINT_CENTROMERE = "point_centromere"
POLYCISTRONIC_PRIMARY_TRANSCRIPT = "polycistronic_primary_transcript"
EXPERIMENTALLY_DEFINED_BINDING_REGION = "experimentally_defined_binding_region"
A_BOX = "A_box"
DIRECT_REPEAT = "direct_repeat"
MATCH_PART = "match_part"
ARCHAEAL_INTRON = "archaeal_intron"
R_THREE_PRIME_LTR_REGION = "R_three_prime_LTR_region"
BETA_TURN = "beta_turn"
DNASEI_HYPERSENSITIVE_SITE = "DNAseI_hypersensitive_site"
METHIONINE_TRNA_PRIMARY_TRANSCRIPT = "methionine_tRNA_primary_transcript"
PROLINE = "proline"
STS = "STS"
INTRAMEMBRANE_POLYPEPTIDE_REGION = "intramembrane_polypeptide_region"
YAC_END = "YAC_end"
V_VDJ_CLUSTER = "V_VDJ_cluster"
OLIGO_U_TAIL = "oligo_U_tail"
JUNCTION = "junction"
INVERSION_BREAKPOINT = "inversion_breakpoint"
MITOCHONDRIAL_TARGETING_SIGNAL = "mitochondrial_targeting_signal"
E_BOX_MOTIF = "E_box_motif"
GENE_SUBARRAY = "gene_subarray"
PROTEIN_PROTEIN_CONTACT = "protein_protein_contact"
PHEROMONE_RESPONSE_ELEMENT = "pheromone_response_element"
TRANSLATIONAL_FRAMESHIFT = "translational_frameshift"
TERMINAL_INVERTED_REPEAT = "terminal_inverted_repeat"
A_BOX_TYPE_1 = "A_box_type_1"
BETA_TURN_RIGHT_HANDED_TYPE_ONE = "beta_turn_right_handed_type_one"
WYBUTOSINE = "wybutosine"
RECODED_MRNA = "recoded_mRNA"
SPLICE_JUNCTION = "splice_junction"
PROPEPTIDE_CLEAVAGE_SITE = "propeptide_cleavage_site"
LIPOPROTEIN_SIGNAL_PEPTIDE = "lipoprotein_signal_peptide"
VIRTUAL_SEQUENCE = "virtual_sequence"
POLYADENYLATED_MRNA = "polyadenylated_mRNA"
ARIA = "ARIA"
AUTOCATALYTICALLY_SPLICED_INTRON = "autocatalytically_spliced_intron"
FLEX_ELEMENT = "FLEX_element"
UPSTREAM_AUG_CODON = "upstream_AUG_codon"
J_GENE_SEGMENT = "J_gene_segment"
AP_1_BINDING_SITE = "AP_1_binding_site"
MORPHOLINO_OLIGO = "morpholino_oligo"
PROMOTER_TRAP_CONSTRUCT = "promoter_trap_construct"
I_MOTIF = "i_motif"
N2_7_DIMETHYLGUANOSINE = "N2_7_dimethylguanosine"
RESTRICTION_ENZYME_BINDING_SITE = "restriction_enzyme_binding_site"
FIVE_PRIME_TERMINAL_INVERTED_REPEAT = "five_prime_terminal_inverted_repeat"
EPITOPE = "epitope"
SRP_RNA = "SRP_RNA"
SERYL_TRNA = "seryl_tRNA"
THREONINE = "threonine"
LIGAND_BINDING_SITE = "ligand_binding_site"
RECOMBINATION_SIGNAL_SEQUENCE = "recombination_signal_sequence"
HISTONE_ACYLATION_REGION = "histone_acylation_region"
NUCLEOMORPHIC_CHROMOSOME = "nucleomorphic_chromosome"
TRNA_REGION = "tRNA_region"
THREE_PRIME_STEM_LOOP_STRUCTURE = "three_prime_stem_loop_structure"
RNA_CHROMOSOME = "RNA_chromosome"
YAC = "YAC"
HISTONE_METHYLATION_SITE = "histone_methylation_site"
INOSINE = "inosine"
GNA_OLIGO = "GNA_oligo"
RECODED_CODON = "recoded_codon"
ANTI_ARRET = "anti_ARRET"
BIOLOGICAL_REGION = "biological_region"
PROMOTER_FLANKING_REGION = "promoter_flanking_region"
GROUP_I_INTRON = "group_I_intron"
NON_CANONICAL_START_CODON = "non_canonical_start_codon"
PARALOGOUS_REGION = "paralogous_region"
TRANSCRIPTION_PAUSE_SITE = "transcription_pause_site"
TERMINATOR_OF_TYPE_2_RNAPOL_III_PROMOTER = "terminator_of_type_2_RNApol_III_promoter"
GENE_WITH_STOP_CODON_READ_THROUGH = "gene_with_stop_codon_read_through"
INSERTION_SITE = "insertion_site"
B_BOX = "B_box"
PAIRED_END_FRAGMENT = "paired_end_fragment"
SERINE_THREONINE_TURN = "serine_threonine_turn"
HISTONE_3_ACETYLATION_SITE = "histone_3_acetylation_site"
NUCLEOTIDE_BINDING_SITE = "nucleotide_binding_site"
N2_7_2PRIRME_O_TRIMETHYLGUANOSINE = "N2_7_2prirme_O_trimethylguanosine"
GAMMA_TURN_INVERSE = "gamma_turn_inverse"
GC_RICH_PROMOTER_REGION = "GC_rich_promoter_region"
DRE = "DRE"
INTERNAL_SHINE_DALGARNO_SEQUENCE = "internal_Shine_Dalgarno_sequence"
AFLP_FRAGMENT = "AFLP_fragment"
SL10_ACCEPTOR_SITE = "SL10_acceptor_site"
TWO_THIOURIDINE = "two_thiouridine"
PYRROLYSINE_TRNA_PRIMARY_TRANSCRIPT = "pyrrolysine_tRNA_primary_transcript"
_5_FORMYLCYTOSINE = "5_formylcytosine"
INTRONIC_LNCRNA = "intronic_lncRNA"
MRNA_READ = "mRNA_read"
NON_CANONICAL_FIVE_PRIME_SPLICE_SITE = "non_canonical_five_prime_splice_site"
H4K16_ACYLATION_SITE = "H4K16_acylation_site"
CODING_EXON = "coding_exon"
INSERTION_SEQUENCE = "insertion_sequence"
PROTEASE_SITE = "protease_site"
RRNA_PRIMARY_TRANSCRIPT_REGION = "rRNA_primary_transcript_region"
MINI_EXON_DONOR_RNA = "mini_exon_donor_RNA"
RECOMBINATION_FEATURE = "recombination_feature"
SERINE = "serine"
REGULON = "regulon"
BASE = "base"
TR_BOX = "TR_box"
BETA_BULGE = "beta_bulge"
VALIDATED_CDNA_CLONE = "validated_cDNA_clone"
J_NONAMER = "J_nonamer"
DELETION_BREAKPOINT = "deletion_breakpoint"
RIBOZYME = "ribozyme"
UTR = "UTR"
GLUTAMINE_TRNA_PRIMARY_TRANSCRIPT = "glutamine_tRNA_primary_transcript"
ONE_TWO_PRIME_O_DIMETHYLADENOSINE = "one_two_prime_O_dimethyladenosine"
ANTISENSE_RNA = "antisense_RNA"
N4_N4_2_PRIME_O_TRIMETHYLCYTIDINE = "N4_N4_2_prime_O_trimethylcytidine"
H_ACA_BOX_SNORNA_PRIMARY_TRANSCRIPT = "H_ACA_box_snoRNA_primary_transcript"
INTERGENIC_REGION = "intergenic_region"
U6ATAC_SNRNA = "U6atac_snRNA"
H3K79_MONOMETHYLATION_SITE = "H3K79_monomethylation_site"
CDS_SUPPORTED_BY_DOMAIN_MATCH_DATA = "CDS_supported_by_domain_match_data"
MIR_ENCODING_TRNA_PRIMARY_TRANSCRIPT = "miR_encoding_tRNA_primary_transcript"
VAULT_RNA = "vault_RNA"
THREE_PRIME_NONCODING_EXON = "three_prime_noncoding_exon"
NESTED_TANDEM_REPEAT = "nested_tandem_repeat"
CLONE = "clone"
R_GNA_OLIGO = "R_GNA_oligo"
METHYLATED_DNA_BASE_FEATURE = "methylated_DNA_base_feature"
ASPARTIC_ACID = "aspartic_acid"
NEGATIVELY_AUTOREGULATED_GENE = "negatively_autoregulated_gene"
TWO_METHYLTHIO_N6_CIS_HYDROXYISOPENTENYL_ADENOSINE = "two_methylthio_N6_cis_hydroxyisopentenyl_adenosine"
C_GENE_SEGMENT = "C_gene_segment"
DART_MARKER = "DArT_marker"
MINOR_TSS = "minor_TSS"
RECOMBINATIONALLY_REARRANGED_VERTEBRATE_IMMUNE_SYSTEM_GENE = "recombinationally_rearranged_vertebrate_immune_system_gene"
H2BK20_ACETYLATION_SITE = "H2BK20_acetylation_site"
C_TO_T_TRANSITION_AT_PCPG_SITE = "C_to_T_transition_at_pCpG_site"
PSEUDOGENE = "pseudogene"
PATHOGENIC_ISLAND = "pathogenic_island"
CHROMOPLAST_CHROMOSOME = "chromoplast_chromosome"
EXTERNAL_TRANSCRIBED_SPACER_REGION = "external_transcribed_spacer_region"
RESTRICTION_ENZYME_REGION = "restriction_enzyme_region"
FIVE_PRIME_D_RECOMBINATION_SIGNAL_SEQUENCE = "five_prime_D_recombination_signal_sequence"
STOP_CODON_READ_THROUGH = "stop_codon_read_through"
ADAPTIVE_ISLAND = "adaptive_island"
H2AK5_ACETYLATION_SITE = "H2AK5_acetylation_site"
MIRTRON = "miRtron"
H3K36_METHYLATION_SITE = "H3K36_methylation_site"
T7_RNA_POLYMERASE_PROMOTER = "T7_RNA_Polymerase_Promoter"
_2A_SELF_CLEAVING_PEPTIDE_REGION = "2A_self_cleaving_peptide_region"
_2KB_DOWNSTREAM_VARIANT = "2KB_downstream_variant"
_2KB_UPSTREAM_VARIANT = "2KB_upstream_variant"
_3D_POLYPEPTIDE_STRUCTURE_VARIANT = "3D_polypeptide_structure_variant"
_3_PRIME_UTR_ELONGATION = "3_prime_UTR_elongation"
_3_PRIME_UTR_EXON_VARIANT = "3_prime_UTR_exon_variant"
_3_PRIME_UTR_INTRON_VARIANT = "3_prime_UTR_intron_variant"
_3_PRIME_UTR_TRUNCATION = "3_prime_UTR_truncation"
_3_PRIME_UTR_VARIANT = "3_prime_UTR_variant"
_500B_DOWNSTREAM_VARIANT = "500B_downstream_variant"
_5KB_DOWNSTREAM_VARIANT = "5KB_downstream_variant"
_5KB_UPSTREAM_VARIANT = "5KB_upstream_variant"
_5S_SINE_RETROTRANSPOSON = "5S_SINE_retrotransposon"
_5_PRIME_UTR_ELONGATION = "5_prime_UTR_elongation"
_5_PRIME_UTR_EXON_VARIANT = "5_prime_UTR_exon_variant"
_5_PRIME_UTR_INTRON_VARIANT = "5_prime_UTR_intron_variant"
_5_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT = "5_prime_UTR_premature_start_codon_gain_variant"
_5_PRIME_UTR_PREMATURE_START_CODON_LOSS_VARIANT = "5_prime_UTR_premature_start_codon_loss_variant"
_5_PRIME_UTR_PREMATURE_START_CODON_VARIANT = "5_prime_UTR_premature_start_codon_variant"
_5_PRIME_UTR_TRUNCATION = "5_prime_UTR_truncation"
_5_PRIME_UTR_UORF_FRAMESHIFT_VARIANT = "5_prime_UTR_uORF_frameshift_variant"
_5_PRIME_UTR_UORF_STOP_CODON_GAIN_VARIANT = "5_prime_UTR_uORF_stop_codon_gain_variant"
_5_PRIME_UTR_UORF_STOP_CODON_LOSS_VARIANT = "5_prime_UTR_uORF_stop_codon_loss_variant"
_5_PRIME_UTR_UORF_STOP_CODON_VARIANT = "5_prime_UTR_uORF_stop_codon_variant"
_5_PRIME_UTR_UORF_VARIANT = "5_prime_UTR_uORF_variant"
_5_PRIME_UTR_VARIANT = "5_prime_UTR_variant"
_7SL_SINE_RETROTRANSPOSON = "7SL_SINE_retrotransposon"
AUG_INITIATED_UORF = "AUG_initiated_uORF"
ALU_DELETION = "Alu_deletion"
ALU_INSERTION = "Alu_insertion"
BEL_PAO_LTR_RETROTRANSPOSON = "Bel_Pao_LTR_retrotransposon"
C_D_H_ACA_BOX_SCARNA = "C-D_H_ACA_box_scaRNA"
C_D_H_ACA_BOX_SCARNA_GENE = "C-D_H_ACA_box_scaRNA_gene"
CAAX_BOX = "CAAX_box"
CACTA_TIR_TRANSPOSON = "CACTA_TIR_transposon"
CARG_BOX = "CArG_box"
CDS_EXTENSION = "CDS_extension"
CDS_FIVE_PRIME_EXTENSION = "CDS_five_prime_extension"
CDS_SUPPORTED_BY_PEPTIDE_SPECTRUM_MATCH = "CDS_supported_by_peptide_spectrum_match"
CDS_THREE_PRIME_EXTENSION = "CDS_three_prime_extension"
C_D_BOX_SCARNA = "C_D_box_scaRNA"
C_D_BOX_SCARNA_GENE = "C_D_box_scaRNA_gene"
C_D_BOX_SNORNA_ENCODING = "C_D_box_snoRNA_encoding"
C_D_BOX_SNORNA_GENE = "C_D_box_snoRNA_gene"
COPIA_LTR_RETROTRANSPOSON = "Copia_LTR_retrotransposon"
CRICK_STRAND = "Crick_strand"
CRYPTON_YR_TRANSPOSON = "Crypton_YR_transposon"
DIRS_YR_RETROTRANSPOSON = "DIRS_YR_retrotransposon"
DNA = "DNA"
DNA_LOOP = "DNA_loop"
DNA_LOOP_ANCHOR = "DNA_loop_anchor"
DNASEI_HYPERSENSITIVE_SITE_CURRENT_SPELLING = "DNaseI_hypersensitive_site"
ENDOGENOUS_RETROVIRUS_LTR_RETROTRANSPOSON = "Endogenous_Retrovirus_LTR_retrotransposon"
FRT_FLANKED = "FRT_flanked"
GNA = "GNA"
GYPSY_LTR_RETROTRANSPOSON = "Gypsy_LTR_retrotransposon"
H2AZK11_ACETYLATION_SITE = "H2AZK11_acetylation_site"
H2AZK13_ACETYLATION_SITE = "H2AZK13_acetylation_site"
H2AZK15_ACETYLATION_SITE = "H2AZK15_acetylation_site"
H2AZK4_ACETYLATION_SITE = "H2AZK4_acetylation_site"
H2AZK7_ACETYLATION_SITE = "H2AZK7_acetylation_site"
H2BK5_ACETYLATION_SITE = "H2BK5_acetylation_site"
H3K23_ACETYLATION_SITE = "H3K23_acetylation_site"
H3K4_ACETYLATION_SITE = "H3K4_acetylation_site"
H4K16_ACETYLATION_SITE = "H4K16_acetylation_site"
H4K5_ACETYLATION_SITE = "H4K5_acetylation_site"
H4K8_ACETYLATION_SITE = "H4K8_acetylation_site"
HERV_DELETION = "HERV_deletion"
HERV_INSERTION = "HERV_insertion"
H_ACA_BOX_SCARNA = "H_ACA_box_scaRNA"
H_ACA_BOX_SCARNA_GENE = "H_ACA_box_scaRNA_gene"
H_ACA_BOX_SNORNA_ENCODING = "H_ACA_box_snoRNA_encoding"
H_ACA_BOX_SNORNA_GENE = "H_ACA_box_snoRNA_gene"
IG_C_GENE = "IG_C_gene"
IG_C_PSEUDOGENE = "IG_C_pseudogene"
IG_D_GENE = "IG_D_gene"
IG_J_GENE = "IG_J_gene"
IG_J_PSEUDOGENE = "IG_J_pseudogene"
IG_V_GENE = "IG_V_gene"
IG_V_PSEUDOGENE = "IG_V_pseudogene"
I_LINE_RETROTRANSPOSON = "I_LINE_retrotransposon"
JOCKEY_LINE_RETROTRANSPOSON = "Jockey_LINE_retrotransposon"
L1_LINE_RETROTRANSPOSON = "L1_LINE_retrotransposon"
LARD = "LARD"
LINE1_DELETION = "LINE1_deletion"
LINE1_INSERTION = "LINE1_insertion"
LNA = "LNA"
LOZ1_RESPONSE_ELEMENT = "LOZ1_response_element"
MNV_ARTIFACT = "MNV_artifact"
MAT2P = "Mat2P"
MAT3M = "Mat3M"
MERLIN_TIR_TRANSPOSON = "Merlin_TIR_transposon"
MUTATOR_TIR_TRANSPOSON = "Mutator_TIR_transposon"
NMD_ESCAPING_VARIANT = "NMD_escaping_variant"
NMD_POLYMORPHIC_PSEUDOGENE_TRANSCRIPT = "NMD_polymorphic_pseudogene_transcript"
NMD_TRANSCRIPT = "NMD_transcript"
NMD_TRANSCRIPT_VARIANT = "NMD_transcript_variant"
NMD_TRIGGERING_VARIANT = "NMD_triggering_variant"
NSD_TRANSCRIPT = "NSD_transcript"
NGARO_YR_RETROTRANSPOSON = "Ngaro_YR_retrotransposon"
PIF_HARBINGER_TIR_TRANSPOSON = "PIF_Harbinger_TIR_transposon"
PNA = "PNA"
P_TIR_TRANSPOSON = "P_TIR_transposon"
PENELOPE_RETROTRANSPOSON = "Penelope_retrotransposon"
PHO7_BINDING_SITE = "Pho7_binding_site"
R2_LINE_RETROTRANSPOSON = "R2_LINE_retrotransposon"
RNA = "RNA"
RNA_7SK = "RNA_7SK"
RNA_7SK_GENE = "RNA_7SK_gene"
RNA_THERMOMETER = "RNA_thermometer"
RTE_LINE_RETROTRANSPOSON = "RTE_LINE_retrotransposon"
R_GNA = "R_GNA"
RETROVIRUS_LTR_RETROTRANSPOSON = "Retrovirus_LTR_retrotransposon"
ROBERTSONIAN_FUSION = "Robertsonian_fusion"
SHP_BOX = "SHP_box"
SNV_ARTIFACT = "SNV_artifact"
SRP_RNA_ENCODING = "SRP_RNA_encoding"
SUMO_INTERACTION_MOTIF = "SUMO_interaction_motif"
SVA_DELETION = "SVA_deletion"
SVA_INSERTION = "SVA_insertion"
S_GNA = "S_GNA"
TFBS_ABLATION = "TFBS_ablation"
TFBS_AMPLIFICATION = "TFBS_amplification"
TFBS_FUSION = "TFBS_fusion"
TFBS_TRANSLOCATION = "TFBS_translocation"
TFRS_COLLECTION = "TFRS_collection"
TFRS_MODULE = "TFRS_module"
TF_BINDING_SITE_VARIANT = "TF_binding_site_variant"
TNA = "TNA"
TRIM = "TRIM"
TR_C_GENE = "TR_C_Gene"
TR_D_GENE = "TR_D_Gene"
TR_J_GENE = "TR_J_Gene"
TR_J_PSEUDOGENE = "TR_J_pseudogene"
TR_V_GENE = "TR_V_Gene"
TR_V_PSEUDOGENE = "TR_V_pseudogene"
T_CELL_RECEPTOR_GENE = "T_cell_receptor_gene"
T_CELL_RECEPTOR_PSEUDOGENE = "T_cell_receptor_pseudogene"
TC1_MARINER_TIR_TRANSPOSON = "Tc1_Mariner_TIR_transposon"
TRANSIB_TIR_TRANSPOSON = "Transib_TIR_transposon"
U14_SNORNA_GENE = "U14_snoRNA_gene"
U3_SNORNA_GENE = "U3_snoRNA_gene"
U7_SNRNA = "U7_snRNA"
UNAAAC_MOTIF = "UNAAAC_motif"
UTR_VARIANT = "UTR_variant"
VIPER_YR_RETROTRANSPOSON = "Viper_YR_retrotransposon"
WATSON_STRAND = "Watson_strand"
YR_RETROTRANSPOSON = "YR_retrotransposon"
Y_RNA_GENE = "Y_RNA_gene"
ZAS1_RECOGNITION_MOTIF = "Zas1_recognition_motif"
ABOLISHED_TRANSCRIPTION_LEVEL_VARIANT = "abolished_transcription_level_variant"
ABSENT_GENE_PRODUCT = "absent_gene_product"
ACCESSIBLE_DNA_REGION = "accessible_DNA_region"
ALLELE = "allele"
ALLELIC_FREQUENCY = "allelic_frequency"
ALLELIC_PSEUDOGENE = "allelic_pseudogene"
ALLELIC_PSEUDOGENIC_RRNA = "allelic_pseudogenic_rRNA"
ALLELIC_PSEUDOGENIC_TRNA = "allelic_pseudogenic_tRNA"
ALLELICALLY_EXCLUDED = "allelically_excluded"
ALLOPOLYPLOID = "allopolyploid"
ALTERATION_ATTRIBUTE = "alteration_attribute"
ALTERED_GENE_PRODUCT_LEVEL = "altered_gene_product_level"
ALTERED_GENE_PRODUCT_SEQUENCE = "altered_gene_product_sequence"
ALTERNATIVELY_SPLICED = "alternatively_spliced"
AMBER_STOP_CODON = "amber_stop_codon"
AMINO_ACID_DELETION = "amino_acid_deletion"
AMINO_ACID_INSERTION = "amino_acid_insertion"
AMINO_ACID_SUBSTITUTION = "amino_acid_substitution"
ANEUPLOID = "aneuploid"
ANEUPLOID_CHROMOSOME = "aneuploid_chromosome"
ANNOTATION_DIRECTED_IMPROVED_DRAFT = "annotation_directed_improved_draft"
ANTISENSE = "antisense"
ANTISENSE_LNCRNA_GENE = "antisense_lncRNA_gene"
APICOPLAST_SEQUENCE = "apicoplast_sequence"
ASSORTMENT_DERIVED_ANEUPLOID = "assortment_derived_aneuploid"
ASSORTMENT_DERIVED_DEFICIENCY = "assortment_derived_deficiency"
ASSORTMENT_DERIVED_DEFICIENCY_PLUS_DUPLICATION = "assortment_derived_deficiency_plus_duplication"
ASSORTMENT_DERIVED_DUPLICATION = "assortment_derived_duplication"
ASSORTMENT_DERIVED_VARIATION = "assortment_derived_variation"
AUTOPOLYPLOID = "autopolyploid"
AUTOREGULATED = "autoregulated"
AUTOSYNAPTIC_CHROMOSOME = "autosynaptic_chromosome"
BACTERIAL_RNAPOL_PROMOTER_SIGMA54_ELEMENT = "bacterial_RNApol_promoter_sigma54_element"
BACTERIAL_RNAPOL_PROMOTER_SIGMA_70_ELEMENT = "bacterial_RNApol_promoter_sigma_70_element"
BACTERIAL_RNAPOL_PROMOTER_SIGMA_ECF_ELEMENT = "bacterial_RNApol_promoter_sigma_ecf_element"
BENIGN_VARIANT = "benign_variant"
BIDIRECTIONAL_GENE_FUSION = "bidirectional_gene_fusion"
BIDIRECTIONAL_PROMOTER_LNCRNA = "bidirectional_promoter_lncRNA"
BIDIRECTIONAL_PROMOTER_LNCRNA_GENE = "bidirectional_promoter_lncRNA_gene"
BIOSYNTHETIC_GENE_CLUSTER = "biosynthetic_gene_cluster"
BIPARTITE_DUPLICATION = "bipartite_duplication"
BIPARTITE_INVERSION = "bipartite_inversion"
BOUND_BY_FACTOR = "bound_by_factor"
BOUND_BY_NUCLEIC_ACID = "bound_by_nucleic_acid"
BOUND_BY_PROTEIN = "bound_by_protein"
CDNA = "cDNA"
CAPPED = "capped"
CASSETTE_ARRAY_MEMBER = "cassette_array_member"
CELLULAR_MISLOCALIZATION_VARIANT = "cellular_mislocalization_variant"
CHLOROPLAST_DNA = "chloroplast_DNA"
CHLOROPLAST_SEQUENCE = "chloroplast_sequence"
CHROMATIN_REGULATORY_REGION = "chromatin_regulatory_region"
CHROMOPLAST_SEQUENCE = "chromoplast_sequence"
CHROMOSOMAL_DELETION = "chromosomal_deletion"
CHROMOSOMAL_DUPLICATION = "chromosomal_duplication"
CHROMOSOMAL_INVERSION = "chromosomal_inversion"
CHROMOSOMAL_TRANSLOCATION = "chromosomal_translocation"
CHROMOSOMAL_TRANSPOSITION = "chromosomal_transposition"
CHROMOSOMAL_VARIATION_ATTRIBUTE = "chromosomal_variation_attribute"
CHROMOSOMALLY_ABERRANT_GENOME = "chromosomally_aberrant_genome"
CHROMOSOME_FISSION = "chromosome_fission"
CHROMOSOME_NUMBER_VARIATION = "chromosome_number_variation"
CHROMOSOME_STRUCTURE_VARIATION = "chromosome_structure_variation"
CHROMOSOME_VARIATION = "chromosome_variation"
CIRCULAR = "circular"
CIRCULAR_MRNA = "circular_mRNA"
CIRCULAR_NCRNA = "circular_ncRNA"
CIRCULAR_PLASMID = "circular_plasmid"
CIS_REGULATORY_MODULE = "cis_regulatory_module"
CODING_SEQUENCE_VARIANT = "coding_sequence_variant"
CODING_TRANSCRIPT_INTRON_VARIANT = "coding_transcript_intron_variant"
CODING_TRANSCRIPT_VARIANT = "coding_transcript_variant"
CODING_TRANSCRIPT_WITH_RETAINED_INTRON = "coding_transcript_with_retained_intron"
CODING_VARIANT_QUALITY = "coding_variant_quality"
CODON_REDEFINED = "codon_redefined"
COMMON_FRAGILE_SITE = "common_fragile_site"
COMMON_VARIANT = "common_variant"
COMPENSATORY_TRANSCRIPT_SECONDARY_STRUCTURE_VARIANT = "compensatory_transcript_secondary_structure_variant"
COMPLEX_3D_STRUCTURAL_VARIANT = "complex_3D_structural_variant"
COMPLEX_CHANGE_OF_TRANSLATIONAL_PRODUCT_VARIANT = "complex_change_of_translational_product_variant"
COMPLEX_CHROMOSOMAL_REARRANGEMENT = "complex_chromosomal_rearrangement"
COMPLEX_OPERON = "complex_operon"
COMPLEX_REGULON = "complex_regulon"
COMPLEX_TRANSCRIPT_VARIANT = "complex_transcript_variant"
COMPOUND_CHROMOSOME = "compound_chromosome"
COMPOUND_CHROMOSOME_ARM = "compound_chromosome_arm"
CONFORMATIONAL_CHANGE_VARIANT = "conformational_change_variant"
CONSENSUS = "consensus"
CONSERVATIVE_AMINO_ACID_SUBSTITUTION = "conservative_amino_acid_substitution"
CONSERVATIVE_INFRAME_DELETION = "conservative_inframe_deletion"
CONSERVATIVE_INFRAME_INSERTION = "conservative_inframe_insertion"
CONSERVATIVE_MISSENSE_VARIANT = "conservative_missense_variant"
CONSERVED = "conserved"
CONSERVED_INTERGENIC_VARIANT = "conserved_intergenic_variant"
CONSERVED_INTRON_VARIANT = "conserved_intron_variant"
COPY_NUMBER_CHANGE = "copy_number_change"
COPY_NUMBER_DECREASE = "copy_number_decrease"
COPY_NUMBER_INCREASE = "copy_number_increase"
CORE_EUKARYOTIC_PROMOTER_ELEMENT = "core_eukaryotic_promoter_element"
CORE_PROKARYOTIC_PROMOTER_ELEMENT = "core_prokaryotic_promoter_element"
CORE_VIRAL_PROMOTER_ELEMENT = "core_viral_promoter_element"
CRYPTIC = "cryptic"
CRYPTIC_PROMOTER = "cryptic_promoter"
CRYPTIC_SPLICE_ACCEPTOR = "cryptic_splice_acceptor"
CRYPTIC_SPLICE_DONOR = "cryptic_splice_donor"
CRYPTIC_SPLICE_SITE_VARIANT = "cryptic_splice_site_variant"
CYANELLE_SEQUENCE = "cyanelle_sequence"
CYCLIC_TRANSLOCATION = "cyclic_translocation"
CYTOSOLIC_16S_RRNA = "cytosolic_16S_rRNA"
CYTOSOLIC_18S_RRNA = "cytosolic_18S_rRNA"
CYTOSOLIC_23S_RRNA = "cytosolic_23S_rRNA"
CYTOSOLIC_25S_RRNA = "cytosolic_25S_rRNA"
CYTOSOLIC_28S_RRNA = "cytosolic_28S_rRNA"
CYTOSOLIC_2S_RRNA = "cytosolic_2S_rRNA"
CYTOSOLIC_5S_RRNA = "cytosolic_5S_rRNA"
CYTOSOLIC_5_8S_RRNA = "cytosolic_5_8S_rRNA"
CYTOSOLIC_LSU_RRNA = "cytosolic_LSU_rRNA"
CYTOSOLIC_LSU_RRNA_GENE = "cytosolic_LSU_rRNA_gene"
CYTOSOLIC_SSU_RRNA = "cytosolic_SSU_rRNA"
CYTOSOLIC_SSU_RRNA_GENE = "cytosolic_SSU_rRNA_gene"
CYTOSOLIC_RRNA = "cytosolic_rRNA"
CYTOSOLIC_RRNA_16S_GENE = "cytosolic_rRNA_16S_gene"
CYTOSOLIC_RRNA_18S_GENE = "cytosolic_rRNA_18S_gene"
CYTOSOLIC_RRNA_23S_GENE = "cytosolic_rRNA_23S_gene"
CYTOSOLIC_RRNA_25S_GENE = "cytosolic_rRNA_25S_gene"
CYTOSOLIC_RRNA_28S_GENE = "cytosolic_rRNA_28S_gene"
CYTOSOLIC_RRNA_2S_GENE = "cytosolic_rRNA_2S_gene"
CYTOSOLIC_RRNA_5S_GENE = "cytosolic_rRNA_5S_gene"
CYTOSOLIC_RRNA_5_8S_GENE = "cytosolic_rRNA_5_8S_gene"
CYTOSOLIC_RRNA_GENE = "cytosolic_rRNA_gene"
DE_NOVO_VARIANT = "de_novo_variant"
DECREASED_ENZYME_ACTIVITY_VARIANT = "decreased_enzyme_activity_variant"
DECREASED_GENE_PRODUCT_LEVEL = "decreased_gene_product_level"
DECREASED_POLYADENYLATION_VARIANT = "decreased_polyadenylation_variant"
DECREASED_POLYPEPTIDE_LEVEL = "decreased_polypeptide_level"
DECREASED_TRANSCRIPT_LEVEL_VARIANT = "decreased_transcript_level_variant"
DECREASED_TRANSCRIPT_STABILITY_VARIANT = "decreased_transcript_stability_variant"
DECREASED_TRANSCRIPTION_RATE_VARIANT = "decreased_transcription_rate_variant"
DECREASED_TRANSLATIONAL_PRODUCT_LEVEL = "decreased_translational_product_level"
DEFICIENT_INTERCHROMOSOMAL_TRANSPOSITION = "deficient_interchromosomal_transposition"
DEFICIENT_INTRACHROMOSOMAL_TRANSPOSITION = "deficient_intrachromosomal_transposition"
DEFICIENT_INVERSION = "deficient_inversion"
DEFICIENT_TRANSLOCATION = "deficient_translocation"
DELETION_ARTIFACT = "deletion_artifact"
DELINS = "delins"
DESIGNED_SEQUENCE = "designed_sequence"
DEXTROSYNAPTIC_CHROMOSOME = "dextrosynaptic_chromosome"
DICISTRONIC = "dicistronic"
DIPLOTYPE = "diplotype"
DIRECT = "direct"
DIRECTION_ATTRIBUTE = "direction_attribute"
DISEASE_ASSOCIATED_VARIANT = "disease_associated_variant"
DISEASE_CAUSING_VARIANT = "disease_causing_variant"
DISRUPTIVE_INFRAME_DELETION = "disruptive_inframe_deletion"
DISRUPTIVE_INFRAME_INSERTION = "disruptive_inframe_insertion"
DOMINANT_NEGATIVE_VARIANT = "dominant_negative_variant"
DOUBLE = "double"
DOUBLE_STRANDED_CDNA = "double_stranded_cDNA"
DOWNSTREAM_GENE_VARIANT = "downstream_gene_variant"
DOWNSTREAM_INTERGENIC_FUSION = "downstream_intergenic_fusion"
DOWNSTREAM_TRANSCRIPT_VARIANT = "downstream_transcript_variant"
DUPLICATION_ARTIFACT = "duplication_artifact"
DUPLICATION_ATTRIBUTE = "duplication_attribute"
EARLY_ORIGIN_OF_REPLICATION = "early_origin_of_replication"
EDITED = "edited"
EDITING_VARIANT = "editing_variant"
ELONGATED_IN_FRAME_POLYPEPTIDE_C_TERMINAL = "elongated_in_frame_polypeptide_C_terminal"
ELONGATED_IN_FRAME_POLYPEPTIDE_N_TERMINAL_ELONGATION = "elongated_in_frame_polypeptide_N_terminal_elongation"
ELONGATED_OUT_OF_FRAME_POLYPEPTIDE_C_TERMINAL = "elongated_out_of_frame_polypeptide_C_terminal"
ELONGATED_OUT_OF_FRAME_POLYPEPTIDE_N_TERMINAL = "elongated_out_of_frame_polypeptide_N_terminal"
ELONGATED_POLYPEPTIDE = "elongated_polypeptide"
ELONGATED_POLYPEPTIDE_C_TERMINAL = "elongated_polypeptide_C_terminal"
ELONGATED_POLYPEPTIDE_N_TERMINAL = "elongated_polypeptide_N_terminal"
ENCODES_1_POLYPEPTIDE = "encodes_1_polypeptide"
ENCODES_ALTERNATE_TRANSCRIPTION_START_SITES = "encodes_alternate_transcription_start_sites"
ENCODES_ALTERNATELY_SPLICED_TRANSCRIPTS = "encodes_alternately_spliced_transcripts"
ENCODES_DIFFERENT_POLYPEPTIDES_DIFFERENT_STOP = "encodes_different_polypeptides_different_stop"
ENCODES_DISJOINT_POLYPEPTIDES = "encodes_disjoint_polypeptides"
ENCODES_GREATER_THAN_1_POLYPEPTIDE = "encodes_greater_than_1_polypeptide"
ENCODES_OVERLAPPING_PEPTIDES = "encodes_overlapping_peptides"
ENCODES_OVERLAPPING_PEPTIDES_DIFFERENT_START = "encodes_overlapping_peptides_different_start"
ENCODES_OVERLAPPING_POLYPEPTIDES_DIFFERENT_START_AND_STOP = "encodes_overlapping_polypeptides_different_start_and_stop"
ENDOGENOUS_RETROVIRAL_SEQUENCE = "endogenous_retroviral_sequence"
ENGINEERED = "engineered"
ENHANCER_BLOCKING_ELEMENT = "enhancer_blocking_element"
ENZYMATIC = "enzymatic"
ENZYMATIC_RNA_GENE = "enzymatic_RNA_gene"
EPIGENETICALLY_MODIFIED = "epigenetically_modified"
EPIGENOMICALLY_MODIFIED_REGION = "epigenomically_modified_region"
EUKARYOTIC_PROMOTER = "eukaryotic_promoter"
EXEMPLAR = "exemplar"
EXON_LOSS_VARIANT = "exon_loss_variant"
EXON_VARIANT = "exon_variant"
EXONIC_SPLICE_REGION_VARIANT = "exonic_splice_region_variant"
EXPERIMENTAL_FEATURE_ATTRIBUTE = "experimental_feature_attribute"
EXPERIMENTALLY_DETERMINED = "experimentally_determined"
EXTENDED_INTRONIC_SPLICE_REGION_VARIANT = "extended_intronic_splice_region_variant"
FEATURE_ABLATION = "feature_ablation"
FEATURE_AMPLIFICATION = "feature_amplification"
FEATURE_ATTRIBUTE = "feature_attribute"
FEATURE_ELONGATION = "feature_elongation"
FEATURE_FUSION = "feature_fusion"
FEATURE_TRANSLOCATION = "feature_translocation"
FEATURE_TRUNCATION = "feature_truncation"
FEATURE_VARIANT = "feature_variant"
FINISHED_GENOME = "finished_genome"
FIVE_PRIME_UTR_PREMATURE_START_CODON_LOCATION_VARIANT = "five_prime_UTR_premature_start_codon_location_variant"
FIVE_PRIME_DUPLICATED_TRANSCRIPT = "five_prime_duplicated_transcript"
FIVE_PRIME_FIVE_PRIME_OVERLAP = "five_prime_five_prime_overlap"
FIVE_PRIME_GENE_FUSION = "five_prime_gene_fusion"
FIVE_PRIME_THREE_PRIME_OVERLAP = "five_prime_three_prime_overlap"
FIXED_VARIANT = "fixed_variant"
FLANKED = "flanked"
FLANKING_REPEAT = "flanking_repeat"
FLOXED = "floxed"
FOREIGN = "foreign"
FORWARD = "forward"
FRAGILE_SITE = "fragile_site"
FRAGMENTARY = "fragmentary"
FRAME_RESTORING_VARIANT = "frame_restoring_variant"
FRAMESHIFT = "frameshift"
FRAMESHIFT_ELONGATION = "frameshift_elongation"
FRAMESHIFT_TRUNCATION = "frameshift_truncation"
FRAMESHIFT_UNIDIRECTIONAL_GENE_FUSION = "frameshift_unidirectional_gene_fusion"
FRAMESHIFT_VARIANT = "frameshift_variant"
FRAMESHIFT_VARIANT_NMD_ESCAPING = "frameshift_variant_NMD_escaping"
FRAMESHIFT_VARIANT_NMD_TRIGGERING = "frameshift_variant_NMD_triggering"
FREE = "free"
FREE_CHROMOSOME_ARM = "free_chromosome_arm"
FREE_DUPLICATION = "free_duplication"
FREE_RING_DUPLICATION = "free_ring_duplication"
FUNCTION_UNCERTAIN_VARIANT = "function_uncertain_variant"
FUNCTIONAL_EFFECT_VARIANT = "functional_effect_variant"
FUNCTIONAL_GENE_REGION = "functional_gene_region"
FUNCTIONALLY_ABNORMAL = "functionally_abnormal"
FUNCTIONALLY_NORMAL = "functionally_normal"
FUSION = "fusion"
GRNA_ENCODING = "gRNA_encoding"
GAIN_OF_FUNCTION_VARIANT = "gain_of_function_variant"
GENE_ARRAY_MEMBER = "gene_array_member"
GENE_ATTRIBUTE = "gene_attribute"
GENE_CASSETTE_MEMBER = "gene_cassette_member"
GENE_FUSION = "gene_fusion"
GENE_SUBARRAY_MEMBER = "gene_subarray_member"
GENE_TO_GENE_FEATURE = "gene_to_gene_feature"
GENE_VARIANT = "gene_variant"
GENIC_DOWNSTREAM_TRANSCRIPT_VARIANT = "genic_downstream_transcript_variant"
GENIC_UPSTREAM_TRANSCRIPT_VARIANT = "genic_upstream_transcript_variant"
GENOME = "genome"
GENOMIC_DNA = "genomic_DNA"
GENOMIC_DNA_CONTIG = "genomic_DNA_contig"
GENOMICALLY_IMPRINTED = "genomically_imprinted"
GENOTYPE = "genotype"
GERMLINE_VARIANT = "germline_variant"
GROUP_IIC_INTRON = "group_IIC_intron"
HAT_TIR_TRANSPOSON = "hAT_TIR_transposon"
HAPLOTYPE = "haplotype"
HETERO_COMPOUND_CHROMOSOME = "hetero_compound_chromosome"
HIGH_QUALITY_DRAFT = "high_quality_draft"
HISTONE_2AZ_ACETYLATION_SITE = "histone_2AZ_acetylation_site"
HISTONE_2A_ACETYLATION_SITE = "histone_2A_acetylation_site"
HISTONE_2B_ACETYLATION_SITE = "histone_2B_acetylation_site"
HISTONE_4_ACETYLATION_SITE = "histone_4_acetylation_site"
HOMO_COMPOUND_CHROMOSOME = "homo_compound_chromosome"
HOMOLOGOUS = "homologous"
HOMOLOGOUS_CHROMOSOME_RECOGNITION_AND_PAIRING_LOCUS = "homologous_chromosome_recognition_and_pairing_locus"
HPRNA = "hpRNA"
HPRNA_GENE = "hpRNA_gene"
HYPERPLOID = "hyperploid"
HYPOPLOID = "hypoploid"
IMMUNOGLOBULIN_GENE = "immunoglobulin_gene"
IMMUNOGLOBULIN_PSEUDOGENE = "immunoglobulin_pseudogene"
IMPRINTING_CONTROL_REGION = "imprinting_control_region"
IMPROVED_HIGH_QUALITY_DRAFT = "improved_high_quality_draft"
INACTIVE_CATALYTIC_SITE = "inactive_catalytic_site"
INACTIVE_LIGAND_BINDING_SITE = "inactive_ligand_binding_site"
INCOMPLETE_TERMINAL_CODON_VARIANT = "incomplete_terminal_codon_variant"
INCOMPLETE_TRANSCRIPT_3UTR_VARIANT = "incomplete_transcript_3UTR_variant"
INCOMPLETE_TRANSCRIPT_5UTR_VARIANT = "incomplete_transcript_5UTR_variant"
INCOMPLETE_TRANSCRIPT_CDS = "incomplete_transcript_CDS"
INCOMPLETE_TRANSCRIPT_CODING_SPLICE_VARIANT = "incomplete_transcript_coding_splice_variant"
INCOMPLETE_TRANSCRIPT_EXONIC_VARIANT = "incomplete_transcript_exonic_variant"
INCOMPLETE_TRANSCRIPT_INTRONIC_VARIANT = "incomplete_transcript_intronic_variant"
INCOMPLETE_TRANSCRIPT_SPLICE_REGION_VARIANT = "incomplete_transcript_splice_region_variant"
INCOMPLETE_TRANSCRIPT_VARIANT = "incomplete_transcript_variant"
INCREASED_ENZYME_ACTIVITY_VARIANT = "increased_enzyme_activity_variant"
INCREASED_GENE_PRODUCT_LEVEL = "increased_gene_product_level"
INCREASED_POLYADENYLATION_VARIANT = "increased_polyadenylation_variant"
INCREASED_TRANSCRIPT_LEVEL_VARIANT = "increased_transcript_level_variant"
INCREASED_TRANSCRIPT_STABILITY_VARIANT = "increased_transcript_stability_variant"
INCREASED_TRANSCRIPTION_RATE_VARIANT = "increased_transcription_rate_variant"
INCREASED_TRANSLATIONAL_PRODUCT_LEVEL = "increased_translational_product_level"
INDEL_ARTIFACT = "indel_artifact"
INDEPENDENTLY_KNOWN = "independently_known"
INERT_DNA_SPACER = "inert_DNA_spacer"
INFRAME = "inframe"
INFRAME_DELETION = "inframe_deletion"
INFRAME_INDEL = "inframe_indel"
INFRAME_INSERTION = "inframe_insertion"
INFRAME_UNIDIRECTIONAL_GENE_FUSION = "inframe_unidirectional_gene_fusion"
INFRAME_VARIANT = "inframe_variant"
INITIATOR_CODON_VARIANT = "initiator_codon_variant"
INSERTION_ARTIFACT = "insertion_artifact"
INSERTION_ATTRIBUTE = "insertion_attribute"
INSERTIONAL = "insertional"
INSERTIONAL_DUPLICATION = "insertional_duplication"
INSIDE_INTRON = "inside_intron"
INSIDE_INTRON_ANTIPARALLEL = "inside_intron_antiparallel"
INSIDE_INTRON_PARALLEL = "inside_intron_parallel"
INTEIN_CONTAINING = "intein_containing"
INTERCHROMOSOMAL = "interchromosomal"
INTERCHROMOSOMAL_DUPLICATION = "interchromosomal_duplication"
INTERCHROMOSOMAL_MUTATION = "interchromosomal_mutation"
INTERCHROMOSOMAL_TRANSLOCATION = "interchromosomal_translocation"
INTERCHROMOSOMAL_TRANSPOSITION = "interchromosomal_transposition"
INTERGENIC_1KB_VARIANT = "intergenic_1kb_variant"
INTERGENIC_VARIANT = "intergenic_variant"
INTERMEDIATE = "intermediate"
INTERNAL_FEATURE_ELONGATION = "internal_feature_elongation"
INTRACHROMOSOMAL = "intrachromosomal"
INTRACHROMOSOMAL_DUPLICATION = "intrachromosomal_duplication"
INTRACHROMOSOMAL_MUTATION = "intrachromosomal_mutation"
INTRACHROMOSOMAL_TRANSLOCATION = "intrachromosomal_translocation"
INTRACHROMOSOMAL_TRANSPOSITION = "intrachromosomal_transposition"
INTRAGENIC_VARIANT = "intragenic_variant"
INTRON_GAIN_VARIANT = "intron_gain_variant"
INTRON_VARIANT = "intron_variant"
INVALIDATED = "invalidated"
INVALIDATED_BY_CHIMERIC_CDNA = "invalidated_by_chimeric_cDNA"
INVALIDATED_BY_GENOMIC_CONTAMINATION = "invalidated_by_genomic_contamination"
INVALIDATED_BY_GENOMIC_POLYA_PRIMED_CDNA = "invalidated_by_genomic_polyA_primed_cDNA"
INVALIDATED_BY_PARTIAL_PROCESSING = "invalidated_by_partial_processing"
INVERSION_ATTRIBUTE = "inversion_attribute"
INVERSION_CUM_TRANSLOCATION = "inversion_cum_translocation"
INVERSION_DERIVED_ANEUPLOID_CHROMOSOME = "inversion_derived_aneuploid_chromosome"
INVERSION_DERIVED_BIPARTITE_DEFICIENCY = "inversion_derived_bipartite_deficiency"
INVERSION_DERIVED_BIPARTITE_DUPLICATION = "inversion_derived_bipartite_duplication"
INVERSION_DERIVED_DEFICIENCY_PLUS_ANEUPLOID = "inversion_derived_deficiency_plus_aneuploid"
INVERSION_DERIVED_DEFICIENCY_PLUS_DUPLICATION = "inversion_derived_deficiency_plus_duplication"
INVERSION_DERIVED_DUPLICATION_PLUS_ANEUPLOID = "inversion_derived_duplication_plus_aneuploid"
INVERTED = "inverted"
INVERTED_INSERTIONAL_DUPLICATION = "inverted_insertional_duplication"
INVERTED_INTERCHROMOSOMAL_TRANSPOSITION = "inverted_interchromosomal_transposition"
INVERTED_INTRACHROMOSOMAL_TRANSPOSITION = "inverted_intrachromosomal_transposition"
INVERTED_RING_CHROMOSOME = "inverted_ring_chromosome"
ISOMIR = "isomiR"
KINETOPLAST = "kinetoplast"
KNOB = "knob"
LAEVOSYNAPTIC_CHROMOSOME = "laevosynaptic_chromosome"
LATE_ORIGIN_OF_REPLICATION = "late_origin_of_replication"
LETHAL_VARIANT = "lethal_variant"
LEUCOPLAST_SEQUENCE = "leucoplast_sequence"
LEVEL_OF_TRANSCRIPT_VARIANT = "level_of_transcript_variant"
LINEAR = "linear"
LINEAR_PLASMID = "linear_plasmid"
LNCRNA_GENE = "lncRNA_gene"
LNCRNA_WITH_RETAINED_INTRON = "lncRNA_with_retained_intron"
LONG_TERMINAL_REPEAT_TRANSCRIPT = "long_terminal_repeat_transcript"
LOOP = "loop"
LOSS_OF_FUNCTION_VARIANT = "loss_of_function_variant"
LOSS_OF_HETEROZYGOSITY = "loss_of_heterozygosity"
LOST_POLYPEPTIDE = "lost_polypeptide"
LOW_COMPLEXITY = "low_complexity"
MRNA_ATTRIBUTE = "mRNA_attribute"
MACRONUCLEAR_SEQUENCE = "macronuclear_sequence"
MATERNAL_VARIANT = "maternal_variant"
MATERNALLY_IMPRINTED = "maternally_imprinted"
MATURE_MIRNA_VARIANT = "mature_miRNA_variant"
MATURE_PROTEIN_REGION_OF_CDS = "mature_protein_region_of_CDS"
MEIOTIC_RECOMBINATION_REGION = "meiotic_recombination_region"
MEMBER_OF_REGULON = "member_of_regulon"
METHYLATION_GUIDE_SNORNA_GENE = "methylation_guide_snoRNA_gene"
MIRNA_ENCODING = "miRNA_encoding"
MICRONUCLEAR_SEQUENCE = "micronuclear_sequence"
MINUS_1_FRAMESHIFT = "minus_1_frameshift"
MINUS_1_FRAMESHIFT_VARIANT = "minus_1_frameshift_variant"
MINUS_1_TRANSLATIONAL_FRAMESHIFT = "minus_1_translational_frameshift"
MINUS_1_TRANSLATIONALLY_FRAMESHIFTED = "minus_1_translationally_frameshifted"
MINUS_2_FRAMESHIFT = "minus_2_frameshift"
MINUS_2_FRAMESHIFT_VARIANT = "minus_2_frameshift_variant"
MINUS_2_TRANSLATIONAL_FRAMESHIFT = "minus_2_translational_frameshift"
MISSENSE_VARIANT = "missense_variant"
MITOCHONDRIAL_DNA = "mitochondrial_DNA"
MITOCHONDRIAL_D_LOOP = "mitochondrial_D_loop"
MITOCHONDRIAL_CONTROL_REGION = "mitochondrial_control_region"
MITOCHONDRIAL_SEQUENCE = "mitochondrial_sequence"
MITOTIC_RECOMBINATION_REGION = "mitotic_recombination_region"
MOBILE = "mobile"
MOBILE_ELEMENT_DELETION = "mobile_element_deletion"
MONOCISTRONIC = "monocistronic"
MORPHOLINO_BACKBONE = "morpholino_backbone"
MT_LSU_RRNA = "mt_LSU_rRNA"
MT_LSU_RRNA_GENE = "mt_LSU_rRNA_gene"
MT_SSU_RRNA = "mt_SSU_rRNA"
MT_SSU_RRNA_GENE = "mt_SSU_rRNA_gene"
MT_RRNA = "mt_rRNA"
MT_RRNA_GENE = "mt_rRNA_gene"
MT_TRNA = "mt_tRNA"
MUTATIONAL_HOTSPOT = "mutational_hotspot"
NATURAL = "natural"
NEGATIVELY_AUTOREGULATED = "negatively_autoregulated"
NO_SEQUENCE_ALTERATION = "no_sequence_alteration"
NON_AUG_INITIATED_UORF = "non_AUG_initiated_uORF"
NON_ALLELIC_HOMOLOGOUS_RECOMBINATION_REGION = "non_allelic_homologous_recombination_region"
NON_CODING_TRANSCRIPT_EXON_VARIANT = "non_coding_transcript_exon_variant"
NON_CODING_TRANSCRIPT_INTRON_VARIANT = "non_coding_transcript_intron_variant"
NON_CODING_TRANSCRIPT_SPLICE_REGION_VARIANT = "non_coding_transcript_splice_region_variant"
NON_CODING_TRANSCRIPT_VARIANT = "non_coding_transcript_variant"
NON_COMPLIMENTARY_STEM = "non_complimentary_stem"
NON_CONSERVATIVE_AMINO_ACID_SUBSTITUTION = "non_conservative_amino_acid_substitution"
NON_CONSERVATIVE_MISSENSE_VARIANT = "non_conservative_missense_variant"
NON_PROTEIN_CODING = "non_protein_coding"
NON_SYNONYMOUS = "non_synonymous"
NONCONTIGUOUS_FINISHED = "noncontiguous_finished"
NONSYNONYMOUS_VARIANT = "nonsynonymous_variant"
NUCLEAR_SEQUENCE = "nuclear_sequence"
NUCLEIC_ACID = "nucleic_acid"
NUCLEOMORPHIC_SEQUENCE = "nucleomorphic_sequence"
NUCLEOTIDE_CLEAVAGE_SITE = "nucleotide_cleavage_site"
NULL_MUTATION = "null_mutation"
OCHRE_STOP_CODON = "ochre_stop_codon"
ONCOGENIC_VARIANT = "oncogenic_variant"
OPAL_STOP_CODON = "opal_stop_codon"
OPERON_MEMBER = "operon_member"
ORGANELLE_SEQUENCE = "organelle_sequence"
ORPHAN = "orphan"
ORTHOLOGOUS = "orthologous"
OVERLAPPING = "overlapping"
PARACENTRIC = "paracentric"
PARACENTRIC_INVERSION = "paracentric_inversion"
PARALOGOUS = "paralogous"
PARTIALLY_CHARACTERIZED_CHROMOSOMAL_MUTATION = "partially_characterized_chromosomal_mutation"
PARTIALLY_DUPLICATED_TRANSCRIPT = "partially_duplicated_transcript"
PATERNAL_VARIANT = "paternal_variant"
PATERNALLY_IMPRINTED = "paternally_imprinted"
PEDIGREE_SPECIFIC_VARIANT = "pedigree_specific_variant"
PEPTIDE_COLLECTION = "peptide_collection"
PEPTIDYL = "peptidyl"
PERICENTRIC = "pericentric"
PERICENTRIC_INVERSION = "pericentric_inversion"
PIGGYBAC_TIR_TRANSPOSON = "piggyBac_TIR_transposon"
PLASMID_LOCATION = "plasmid_location"
PLASTID_LSU_RRNA = "plastid_LSU_rRNA"
PLASTID_LSU_RRNA_GENE = "plastid_LSU_rRNA_gene"
PLASTID_SSU_RRNA = "plastid_SSU_rRNA"
PLASTID_SSU_RRNA_GENE = "plastid_SSU_rRNA_gene"
PLASTID_RRNA = "plastid_rRNA"
PLASTID_RRNA_GENE = "plastid_rRNA_gene"
PLASTID_SEQUENCE = "plastid_sequence"
PLUS_1_FRAMESHIFT = "plus_1_frameshift"
PLUS_1_FRAMESHIFT_VARIANT = "plus_1_frameshift_variant"
PLUS_1_TRANSLATIONALLY_FRAMESHIFTED = "plus_1_translationally_frameshifted"
PLUS_2_FRAMESHIFT_VARIANT = "plus_2_frameshift_variant"
PLUS_2_FRAMSHIFT = "plus_2_framshift"
POLYA_SITE_CLUSTER = "polyA_site_cluster"
POLYADENYLATED = "polyadenylated"
POLYADENYLATION_VARIANT = "polyadenylation_variant"
POLYCISTRONIC = "polycistronic"
POLYMER_ATTRIBUTE = "polymer_attribute"
POLYMORPHIC_PSEUDOGENE_PROCESSED_TRANSCRIPT = "polymorphic_pseudogene_processed_transcript"
POLYMORPHIC_PSEUDOGENE_WITH_RETAINED_INTRON = "polymorphic_pseudogene_with_retained_intron"
POLYMORPHIC_SEQUENCE_VARIANT = "polymorphic_sequence_variant"
POLYMORPHIC_VARIANT = "polymorphic_variant"
POLYPEPTIDE_FUNCTION_VARIANT = "polypeptide_function_variant"
POLYPEPTIDE_FUSION = "polypeptide_fusion"
POLYPEPTIDE_GAIN_OF_FUNCTION_VARIANT = "polypeptide_gain_of_function_variant"
POLYPEPTIDE_LOCALIZATION_VARIANT = "polypeptide_localization_variant"
POLYPEPTIDE_LOSS_OF_FUNCTION_VARIANT = "polypeptide_loss_of_function_variant"
POLYPEPTIDE_PARTIAL_LOSS_OF_FUNCTION = "polypeptide_partial_loss_of_function"
POLYPEPTIDE_POST_TRANSLATIONAL_PROCESSING_VARIANT = "polypeptide_post_translational_processing_variant"
POLYPEPTIDE_SEQUENCE_VARIANT = "polypeptide_sequence_variant"
POLYPEPTIDE_TRUNCATION = "polypeptide_truncation"
POLYPLOID = "polyploid"
POPULATION_SPECIFIC_VARIANT = "population_specific_variant"
POSITIVELY_AUTOREGULATED = "positively_autoregulated"
POST_TRANSLATIONALLY_REGULATED = "post_translationally_regulated"
POST_TRANSLATIONALLY_REGULATED_BY_PROTEIN_MODIFICATION = "post_translationally_regulated_by_protein_modification"
POST_TRANSLATIONALLY_REGULATED_BY_PROTEIN_STABILITY = "post_translationally_regulated_by_protein_stability"
PRECURSOR_RNA = "precursor_RNA"
PREDICTED = "predicted"
PREDICTED_BY_AB_INITIO_COMPUTATION = "predicted_by_ab_initio_computation"
PREDICTED_TRANSCRIPT = "predicted_transcript"
PRESENCE_ABSENCE_VARIATION = "presence_absence_variation"
PROCESSED_PSEUDOGENIC_RRNA = "processed_pseudogenic_rRNA"
PROCESSED_PSEUDOGENIC_TRNA = "processed_pseudogenic_tRNA"
PROKARYOTIC_PROMOTER = "prokaryotic_promoter"
PROPEPTIDE_REGION_OF_CDS = "propeptide_region_of_CDS"
PROPLASTID_SEQUENCE = "proplastid_sequence"
PROTEIN_ALTERING_VARIANT = "protein_altering_variant"
PROTEIN_CODING = "protein_coding"
PROVIRAL_LOCATION = "proviral_location"
PSEUDOEXON_VARIANT = "pseudoexon_variant"
PSEUDOGENE_PROCESSED_TRANSCRIPT = "pseudogene_processed_transcript"
PSEUDOGENIC_CDS = "pseudogenic_CDS"
PSEUDOGENIC_TRANSCRIPT_WITH_RETAINED_INTRON = "pseudogenic_transcript_with_retained_intron"
PSEUDOURIDYLATION_GUIDE_SNORNA_GENE = "pseudouridylation_guide_snoRNA_gene"
PUMILIO_RESPONSE_ELEMENT = "pumilio_response_element"
PYRROLYSINE_LOSS = "pyrrolysine_loss"
QUALITY_VALUE = "quality_value"
QUANTITATIVE_VARIANT = "quantitative_variant"
RRNA_ENCODING = "rRNA_encoding"
RANDOM_SEQUENCE = "random_sequence"
RANGE_EXTENDER_ELEMENT = "range_extender_element"
RARE_AMINO_ACID_VARIANT = "rare_amino_acid_variant"
RARE_FRAGILE_SITE = "rare_fragile_site"
RARE_VARIANT = "rare_variant"
RATE_OF_TRANSCRIPTION_VARIANT = "rate_of_transcription_variant"
REARRANGED_AT_DNA_LEVEL = "rearranged_at_DNA_level"
REARRANGEMENT_REGION = "rearrangement_region"
RECIPROCAL = "reciprocal"
RECIPROCAL_CHROMOSOMAL_TRANSLOCATION = "reciprocal_chromosomal_translocation"
RECODED = "recoded"
RECODED_BY_TRANSLATIONAL_BYPASS = "recoded_by_translational_bypass"
RECOMBINATIONALLY_REARRANGED = "recombinationally_rearranged"
REDUNDANT_INSERTED_STOP_GAINED = "redundant_inserted_stop_gained"
REF_MIRNA = "ref_miRNA"
REFERENCE_GENOME = "reference_genome"
REGULATED = "regulated"
REGULATORY_REGION_ABLATION = "regulatory_region_ablation"
REGULATORY_REGION_AMPLIFICATION = "regulatory_region_amplification"
REGULATORY_REGION_FUSION = "regulatory_region_fusion"
REGULATORY_REGION_TRANSLOCATION = "regulatory_region_translocation"
REGULATORY_REGION_VARIANT = "regulatory_region_variant"
REPEAT_INSTABILITY_REGION = "repeat_instability_region"
REPLICATION_START_SITE = "replication_start_site"
RESCUE = "rescue"
RESPONSE_ELEMENT = "response_element"
RETROTRANSPOSED = "retrotransposed"
REVERSE = "reverse"
RIBOZYME_GENE = "ribozyme_gene"
RIBOZYMIC = "ribozymic"
RING_CHROMOSOME = "ring_chromosome"
SBRNA = "sbRNA"
SBRNA_GENE = "sbRNA_gene"
SCRNA_ENCODING = "scRNA_encoding"
SCARNA = "scaRNA"
SCARNA_GENE = "scaRNA_gene"
SCORE = "score"
SELECTION_MARKER = "selection_marker"
SELENOCYSTEINE_GAIN = "selenocysteine_gain"
SELENOCYSTEINE_LOSS = "selenocysteine_loss"
SELF_CLEAVING_RIBOZYME = "self_cleaving_ribozyme"
SENSE_INTRONIC_LNCRNA = "sense_intronic_lncRNA"
SENSE_INTRONIC_LNCRNA_GENE = "sense_intronic_lncRNA_gene"
SENSE_OVERLAP_LNCRNA = "sense_overlap_lncRNA"
SENSE_OVERLAP_LNCRNA_GENE = "sense_overlap_lncRNA_gene"
SEQUENCE_ALTERATION_ARTIFACT = "sequence_alteration_artifact"
SEQUENCE_ATTRIBUTE = "sequence_attribute"
SEQUENCE_COLLECTION = "sequence_collection"
SEQUENCE_COMPARISON = "sequence_comparison"
SEQUENCE_LENGTH_ALTERATION = "sequence_length_alteration"
SEQUENCE_LENGTH_VARIANT = "sequence_length_variant"
SEQUENCE_LOCATION = "sequence_location"
SEQUENCE_SOURCE = "sequence_source"
SEQUENCE_VARIANT = "sequence_variant"
SHORT_TANDEM_REPEAT_CHANGE = "short_tandem_repeat_change"
SHORT_TANDEM_REPEAT_CONTRACTION = "short_tandem_repeat_contraction"
SHORT_TANDEM_REPEAT_EXPANSION = "short_tandem_repeat_expansion"
SHORT_TANDEM_REPEAT_VARIATION = "short_tandem_repeat_variation"
SIGNAL_PEPTIDE_REGION_OF_CDS = "signal_peptide_region_of_CDS"
SILENCED = "silenced"
SILENCED_BY_DNA_METHYLATION = "silenced_by_DNA_methylation"
SILENCED_BY_DNA_MODIFICATION = "silenced_by_DNA_modification"
SILENCED_BY_RNA_INTERFERENCE = "silenced_by_RNA_interference"
SILENCED_BY_HISTONE_DEACETYLATION = "silenced_by_histone_deacetylation"
SILENCED_BY_HISTONE_METHYLATION = "silenced_by_histone_methylation"
SILENCED_BY_HISTONE_MODIFICATION = "silenced_by_histone_modification"
SILENT_MUTATION = "silent_mutation"
SIMPLE_OPERON = "simple_operon"
SIMPLE_REGULON = "simple_regulon"
SINGLE = "single"
SINGLE_STRANDED_CDNA = "single_stranded_cDNA"
SISRNA = "sisRNA"
SNRNA_ENCODING = "snRNA_encoding"
SNCRNA = "sncRNA"
SNCRNA_GENE = "sncRNA_gene"
SNORNA_ENCODING = "snoRNA_encoding"
SOMATIC_VARIANT = "somatic_variant"
SPLICE_ACCEPTOR_VARIANT = "splice_acceptor_variant"
SPLICE_ACCEPTOR_VARIANT_NMD_ESCAPING = "splice_acceptor_variant_NMD_escaping"
SPLICE_ACCEPTOR_VARIANT_NMD_TRIGGERING = "splice_acceptor_variant_NMD_triggering"
SPLICE_DONOR_5TH_BASE_VARIANT = "splice_donor_5th_base_variant"
SPLICE_DONOR_REGION_VARIANT = "splice_donor_region_variant"
SPLICE_DONOR_VARIANT = "splice_donor_variant"
SPLICE_DONOR_VARIANT_NMD_ESCAPING = "splice_donor_variant_NMD_escaping"
SPLICE_DONOR_VARIANT_NMD_TRIGGERING = "splice_donor_variant_NMD_triggering"
SPLICE_POLYPYRIMIDINE_TRACT_VARIANT = "splice_polypyrimidine_tract_variant"
SPLICE_REGION_VARIANT = "splice_region_variant"
SPLICE_SITE_VARIANT = "splice_site_variant"
SPLICING_VARIANT = "splicing_variant"
SPURIOUS_PROTEIN = "spurious_protein"
STRNA_ENCODING = "stRNA_encoding"
STANDARD_DRAFT = "standard_draft"
START_LOST = "start_lost"
START_RETAINED_VARIANT = "start_retained_variant"
STATUS = "status"
STEM = "stem"
STEM_LOOP_REGION = "stem_loop_region"
STOP_GAINED = "stop_gained"
STOP_GAINED_NMD_ESCAPING = "stop_gained_NMD_escaping"
STOP_GAINED_NMD_TRIGGERING = "stop_gained_NMD_triggering"
STOP_LOST = "stop_lost"
STOP_RETAINED_VARIANT = "stop_retained_variant"
STRAND_ATTRIBUTE = "strand_attribute"
STRUCTURAL_INTERACTION_VARIANT = "structural_interaction_variant"
STRUCTURAL_VARIANT = "structural_variant"
SUBSTITUTION_ARTIFACT = "substitution_artifact"
SUPPORTED_BY_EST_OR_CDNA = "supported_by_EST_or_cDNA"
SUPPORTED_BY_DOMAIN_MATCH = "supported_by_domain_match"
SUPPORTED_BY_SEQUENCE_SIMILARITY = "supported_by_sequence_similarity"
SYNONYMOUS = "synonymous"
SYNONYMOUS_VARIANT = "synonymous_variant"
SYNTENIC = "syntenic"
SYNTHETIC_SEQUENCE = "synthetic_sequence"
TRNA_SINE_RETROTRANSPOSON = "tRNA_SINE_retrotransposon"
TRNA_ENCODING = "tRNA_encoding"
TANDEM = "tandem"
TEB1_RECOGNITION_MOTIF = "teb1_recognition_motif"
TELOMERIC_D_LOOP = "telomeric_D_loop"
TERMINAL_REPEAT = "terminal_repeat"
TERMINATOR_CODON_VARIANT = "terminator_codon_variant"
THREE_PRIME_DUPLICATED_TRANSCRIPT = "three_prime_duplicated_transcript"
THREE_PRIME_FIVE_PRIME_OVERLAP = "three_prime_five_prime_overlap"
THREE_PRIME_GENE_FUSION = "three_prime_gene_fusion"
THREE_PRIME_OVERLAPPING_NCRNA = "three_prime_overlapping_ncrna"
THREE_PRIME_THREE_PRIME_OVERLAP = "three_prime_three_prime_overlap"
TMRNA_ENCODING = "tmRNA_encoding"
TOPOLOGICALLY_ASSOCIATED_DOMAIN = "topologically_associated_domain"
TOPOLOGICALLY_ASSOCIATED_DOMAIN_BOUNDARY = "topologically_associated_domain_boundary"
TOPOLOGY_ATTRIBUTE = "topology_attribute"
TRANS_SPLICED = "trans_spliced"
TRANSCRIBED_PROCESSED_PSEUDOGENE = "transcribed_processed_pseudogene"
TRANSCRIBED_UNITARY_PSEUDOGENE = "transcribed_unitary_pseudogene"
TRANSCRIBED_UNPROCESSED_PSEUDOGENE = "transcribed_unprocessed_pseudogene"
TRANSCRIPT_ABLATION = "transcript_ablation"
TRANSCRIPT_AMPLIFICATION = "transcript_amplification"
TRANSCRIPT_ATTRIBUTE = "transcript_attribute"
TRANSCRIPT_FUNCTION_VARIANT = "transcript_function_variant"
TRANSCRIPT_FUSION = "transcript_fusion"
TRANSCRIPT_PROCESSING_VARIANT = "transcript_processing_variant"
TRANSCRIPT_REGULATORY_REGION_FUSION = "transcript_regulatory_region_fusion"
TRANSCRIPT_SECONDARY_STRUCTURE_VARIANT = "transcript_secondary_structure_variant"
TRANSCRIPT_STABILITY_VARIANT = "transcript_stability_variant"
TRANSCRIPT_TRANSLOCATION = "transcript_translocation"
TRANSCRIPT_VARIANT = "transcript_variant"
TRANSCRIPTION_FACTOR_REGULATORY_SITE = "transcription_factor_regulatory_site"
TRANSCRIPTION_TERMINATION_SIGNAL = "transcription_termination_signal"
TRANSCRIPTION_UNIT = "transcription_unit"
TRANSCRIPTION_VARIANT = "transcription_variant"
TRANSCRIPTIONALLY_CONSTITUTIVE = "transcriptionally_constitutive"
TRANSCRIPTIONALLY_INDUCED = "transcriptionally_induced"
TRANSCRIPTIONALLY_REGULATED = "transcriptionally_regulated"
TRANSCRIPTIONALLY_REPRESSED = "transcriptionally_repressed"
TRANSGENIC = "transgenic"
TRANSIT_PEPTIDE_REGION_OF_CDS = "transit_peptide_region_of_CDS"
TRANSLATED_PROCESSED_PSEUDOGENE = "translated_processed_pseudogene"
TRANSLATED_UNPROCESSED_PSEUDOGENE = "translated_unprocessed_pseudogene"
TRANSLATIONAL_PRODUCT_FUNCTION_VARIANT = "translational_product_function_variant"
TRANSLATIONAL_PRODUCT_LEVEL_VARIANT = "translational_product_level_variant"
TRANSLATIONAL_PRODUCT_STRUCTURE_VARIANT = "translational_product_structure_variant"
TRANSLATIONALLY_FRAMESHIFTED = "translationally_frameshifted"
TRANSLATIONALLY_REGULATED = "translationally_regulated"
TRANSLOCATION_ELEMENT = "translocation_element"
TRANSLOCATON_ATTRIBUTE = "translocaton_attribute"
TRINUCLEOTIDE_REPEAT_EXPANSION = "trinucleotide_repeat_expansion"
UNCHARACTERIZED_CHROMOSOMAL_MUTATION = "uncharacterized_chromosomal_mutation"
UNCONFIRMED_TRANSCRIPT = "unconfirmed_transcript"
UNIDIRECTIONAL_GENE_FUSION = "unidirectional_gene_fusion"
UNIDIRECTIONAL_INTERGENIC_FUSION = "unidirectional_intergenic_fusion"
UNINVERTED_INSERTIONAL_DUPLICATION = "uninverted_insertional_duplication"
UNINVERTED_INTERCHROMOSOMAL_TRANSPOSITION = "uninverted_interchromosomal_transposition"
UNINVERTED_INTRACHROMOSOMAL_TRANSPOSITION = "uninverted_intrachromosomal_transposition"
UNIQUE_VARIANT = "unique_variant"
UNIT_OF_GENE_EXPRESSION = "unit_of_gene_expression"
UNITARY_PSEUDOGENIC_RRNA = "unitary_pseudogenic_rRNA"
UNITARY_PSEUDOGENIC_TRNA = "unitary_pseudogenic_tRNA"
UNORIENTED_INSERTIONAL_DUPLICATION = "unoriented_insertional_duplication"
UNORIENTED_INTERCHROMOSOMAL_TRANSPOSITION = "unoriented_interchromosomal_transposition"
UNORIENTED_INTRACHROMOSOMAL_TRANSPOSITION = "unoriented_intrachromosomal_transposition"
UNPROCESSED_PSEUDOGENIC_RRNA = "unprocessed_pseudogenic_rRNA"
UNPROCESSED_PSEUDOGENIC_TRNA = "unprocessed_pseudogenic_tRNA"
UNSPECIFIED_INDEL = "unspecified_indel"
UPSTREAM_GENE_VARIANT = "upstream_gene_variant"
UPSTREAM_INTERGENIC_FUSION = "upstream_intergenic_fusion"
UPSTREAM_TRANSCRIPT_VARIANT = "upstream_transcript_variant"
VALIDATED = "validated"
VARIANT_COLLECTION = "variant_collection"
VARIANT_FREQUENCY = "variant_frequency"
VARIANT_GENOME = "variant_genome"
VARIANT_ORIGIN = "variant_origin"
VARIANT_PHENOTYPE = "variant_phenotype"
VARIANT_QUALITY = "variant_quality"
VAULT_RNA_GENE = "vault_RNA_gene"
VERTEBRATE_IMMUNE_SYSTEM_GENE = "vertebrate_immune_system_gene"
VERTEBRATE_IMMUNE_SYSTEM_PSEUDOGENE = "vertebrate_immune_system_pseudogene"
VIRAL_PROMOTER = "viral_promoter"
WHOLE_GENOME_SEQUENCE_STATUS = "whole_genome_sequence_status"
WILD_TYPE = "wild_type"

# every SO term; a set rather than an enum, being only used to look up GFF3 types
SO_SEQUENCE_TERMS = frozenset({HISTIDYL_TRNA,
                               C_D_BOX_SNORNA_PRIMARY_TRANSCRIPT,
                               POLYCISTRONIC_MRNA,
                               RRNA_23S,
                               TYROSINE_TRNA_PRIMARY_TRANSCRIPT,
                               SHRNA_PRIMARY_TRANSCRIPT,
                               CONTIG,
                               U3_SNORNA,
                               SUPERCONTIG,
                               CIRCULAR_SINGLE_STRANDED_DNA_CHROMOSOME,
                               RNA_APTAMER,
                               ASX_TURN_RIGHT_HANDED_TYPE_ONE,
                               RECOMBINATION_HOTSPOT,
                               GENE_SEGMENT,
                               J_GENE_RECOMBINATION_FEATURE,
                               POLYPEPTIDE_LIGAND_CONTACT,
                               OXYS_RNA,
                               MRNA_WITH_PLUS_2_FRAMESHIFT,
                               LOW_COMPLEXITY_REGION,
                               TARGETING_VECTOR,
                               GENETIC_MARKER,
                               CONSENSUS_AFLP_FRAGMENT,
                               BREU_MOTIF,
                               GENE_WITH_TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT,
                               COPY_NUMBER_VARIATION,
                               EDITED_TRANSCRIPT,
                               INTERNAL_ELIMINATED_SEQUENCE,
                               ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT,
                               LYSOSOMAL_LOCALIZATION_SIGNAL,
                               FIVE_PRIME_UST,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA_70,
                               VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_REARRANGED_SEGMENT,
                               DSR_MOTIF,
                               FIVE_PRIME_INTRON,
                               GENE_TRAP_CONSTRUCT,
                               LOR,
                               ABERRANT_PROCESSED_TRANSCRIPT,
                               DMV4_MOTIF,
                               V_NONAMER,
                               STEROL_REGULATORY_ELEMENT,
                               FIVE_PRIME_D_HEPTAMER,
                               PEROXYWYBUTOSINE,
                               PRIMARY_TRANSCRIPT_REGION,
                               LAMBDA_VECTOR,
                               TWO_PRIME_O_METHYLPSEUDOURIDINE,
                               INTERNAL_UTR,
                               TILING_PATH_CLONE,
                               BACTERIAL_RNAPOL_PROMOTER,
                               DS_OLIGO,
                               TARGET_SITE_DUPLICATION,
                               DPE_MOTIF,
                               SHADOW_ENHANCER,
                               CRM,
                               NEGATIVE_SENSE_SSRNA_VIRAL_SEQUENCE,
                               HOMOL_E_BOX,
                               D_DJ_CLUSTER,
                               FIVE_PRIME_D_SPACER,
                               THREE_TWO_PRIME_O_DIMETHYLURIDINE,
                               EST_MATCH,
                               GLUTAMINE,
                               METHIONYL_TRNA,
                               PHENYLALANINE,
                               TRANSMEMBRANE_POLYPEPTIDE_REGION,
                               CHROMOSOMAL_REGULATORY_ELEMENT,
                               ARGININE,
                               ENDONUCLEASE_SPLICED_INTRON,
                               TRANSPOSABLE_ELEMENT_INSERTION_SITE,
                               MACRONUCLEUS_DESTINED_SEGMENT,
                               V_VJ_C_CLUSTER,
                               MODIFIED_CYTOSINE,
                               CODING_START,
                               TWO_METHYLTHIO_N6_HYDROXYNORVALYL_CARBAMOYLADENOSINE,
                               U3_THREE_PRIME_LTR_REGION,
                               VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_GENE_CLUSTER,
                               CIS_ACTING_HOMOLOGOUS_CHROMOSOME_PAIRING_REGION,
                               U3_LTR_REGION,
                               READING_FRAME,
                               SUBTELOMERE,
                               INVALIDATED_CDNA_CLONE,
                               SPLICED_LEADER_RNA,
                               TILING_PATH_FRAGMENT,
                               IRLINV_SITE,
                               PROVIRAL_REGION,
                               TRANSPOSABLE_ELEMENT_FLANKING_REGION,
                               TRANSCRIPTION_START_CLUSTER,
                               V_GENE_SEGMENT,
                               EXON_JUNCTION,
                               SAP1_RECOGNITION_MOTIF,
                               DISABLED_READING_FRAME,
                               DEFECTIVE_CONJUGATIVE_TRANSPOSON,
                               T_TO_A_TRANSVERSION,
                               VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_SIGNAL_FEATURE,
                               THREE_PRIME_UTR,
                               STRUCTURAL_ALTERATION,
                               CAPPED_MRNA,
                               PROTEIN_CODING_GENE,
                               H3K36_MONOMETHYLATION_SITE,
                               _8_OXOADENINE,
                               SPLICING_REGULATORY_REGION,
                               HOMOL_D_BOX,
                               CENTROMERE_DNA_ELEMENT_I,
                               EST,
                               ENGINEERED_FOREIGN_GENE,
                               EXPRESSED_SEQUENCE_MATCH,
                               TRANSCRIPT_BOUND_BY_PROTEIN,
                               ATTB_SITE,
                               PRE_EDITED_MRNA,
                               DICISTRONIC_PRIMARY_TRANSCRIPT,
                               MODIFIED_L_METHIONINE,
                               UAA_STOP_CODON_SIGNAL,
                               NATURAL_TRANSPOSABLE_ELEMENT,
                               FIVE_PRIME_LTR_COMPONENT,
                               PAC,
                               P_ELEMENT,
                               CENTROMERE_DNA_ELEMENT_III,
                               LINCRNA,
                               BETA_BULGE_LOOP,
                               RPRA_RNA,
                               ANCHOR_BINDING_SITE,
                               D_DJ_J_C_CLUSTER,
                               TWO_PRIME_O_METHYLURIDINE,
                               Y_RNA_PRIMARY_TRANSCRIPT,
                               A_MINOR_RNA_MOTIF,
                               AU_RICH_ELEMENT,
                               FIVE_PRIME_OPEN_READING_FRAME,
                               UNIGENE_CLUSTER,
                               POLYPEPTIDE_REPEAT,
                               RNA_JUNCTION_LOOP,
                               COINTEGRATED_PLASMID,
                               CLASS_II_RNA,
                               IRRINV_SITE,
                               RRNA_LARGE_SUBUNIT_PRIMARY_TRANSCRIPT,
                               CDNA_MATCH,
                               NDM3_MOTIF,
                               CLONE_INSERT_END,
                               CENTROMERE,
                               SELENOCYSTEINYL_TRNA,
                               DG_REPEAT,
                               PROLINE_TRNA_PRIMARY_TRANSCRIPT,
                               TRANSCRIBED_SPACER_REGION,
                               INTEGRATED_PLASMID,
                               CATMAT_RIGHT_HANDED_THREE,
                               REGIONAL_CENTROMERE_OUTER_REPEAT_TRANSCRIPT,
                               SEVEN_METHYLGUANOSINE,
                               C_TO_T_TRANSITION,
                               POLYPEPTIDE_DNA_CONTACT,
                               NON_TERMINAL_RESIDUE,
                               SEQUENCE_DIFFERENCE,
                               REPEAT_FRAGMENT,
                               DMV3_MOTIF,
                               TSS_REGION,
                               RESTRICTION_ENZYME_SINGLE_STRAND_OVERHANG,
                               ENHANCER_TRAP_CONSTRUCT,
                               PSEUDOGENE_BY_UNEQUAL_CROSSING_OVER,
                               ALPHA_BETA_MOTIF,
                               APICOPLAST_CHROMOSOME,
                               N6_ACETYLADENOSINE,
                               MODIFIED_L_GLUTAMIC_ACID,
                               FIVE_PRIME_CIS_SPLICE_SITE,
                               SYNTHETIC_OLIGO,
                               NON_ADJACENT_RESIDUES,
                               CAP,
                               ASSEMBLY,
                               CYSTEINE,
                               GROUP_II_INTRON,
                               RECOMBINATION_REGULATORY_REGION,
                               MODIFIED_L_ASPARAGINE,
                               SECIS_ELEMENT,
                               U5_THREE_PRIME_LTR_REGION,
                               REVERSE_HOOGSTEEN_BASE_PAIR,
                               OLIGO,
                               SMFISH_PROBE,
                               H_ACA_BOX_SNORNA,
                               V_D_DJ_CLUSTER,
                               PARTIAL_GENOMIC_SEQUENCE_ASSEMBLY,
                               PSEUDOGENIC_REGION,
                               CONSENSUS_MRNA,
                               HELITRON,
                               SEVEN_DEAZAGUANOSINE,
                               THREE_PRIME_LTR,
                               CHROMOSOME_BREAKPOINT,
                               SNRNA,
                               SNORNA_GENE,
                               H4K_ACYLATION_REGION,
                               I_BOX,
                               H3K27_TRIMETHYLATION_SITE,
                               TRANS_SPLICE_ACCEPTOR_SITE,
                               LINEAR_SINGLE_STRANDED_RNA_CHROMOSOME,
                               ASPE_PRIMER,
                               PNA_OLIGO,
                               MITE,
                               SEQUENCE_MOTIF,
                               ASPARAGINYL_TRNA,
                               THREONYL_TRNA,
                               TYROSYL_TRNA,
                               SELENOCYSTEINE,
                               MRNA_WITH_MINUS_2_FRAMESHIFT,
                               GENE_GROUP_REGULATORY_REGION,
                               ENHANCER_BINDING_SITE,
                               MRNA_RECODED_BY_TRANSLATIONAL_BYPASS,
                               BASE_CALL_ERROR_CORRECTION,
                               CDS_INDEPENDENTLY_KNOWN,
                               BETA_TURN_LEFT_HANDED_TYPE_TWO,
                               GAMMA_TURN,
                               RESTRICTION_FRAGMENT,
                               SP6_RNA_POLYMERASE_PROMOTER,
                               TRANSCRIPT_REGION,
                               SYMMETRIC_RNA_INTERNAL_LOOP,
                               MINUS_24_SIGNAL,
                               ONE_TWO_PRIME_O_DIMETHYLINOSINE,
                               POLYPEPTIDE_MANGANESE_ION_CONTACT_SITE,
                               CHIMERIC_CDNA_CLONE,
                               PARTIALLY_PROCESSED_CDNA_CLONE,
                               INVERSION,
                               HEPTAMER_OF_RECOMBINATION_FEATURE_OF_VERTEBRATE_IMMUNE_SYSTEM_GENE,
                               FIVE_HYDROXYURIDINE,
                               V_J_C_CLUSTER,
                               H3K9_DIMETHYLATION_SITE,
                               CSL_RESPONSE_ELEMENT,
                               ST_TURN_LEFT_HANDED_TYPE_TWO,
                               CRE,
                               POLYPEPTIDE_MOLYBDENUM_ION_CONTACT_SITE,
                               G_BOX,
                               TRANS_SPLICE_JUNCTION,
                               HELIX_TURN_HELIX,
                               GENE_GROUP,
                               PRIMER_MATCH,
                               PRIRNA,
                               FIVE_PRIME_UTR_INTRON,
                               CDS_SUPPORTED_BY_SEQUENCE_SIMILARITY_DATA,
                               CHROMOSOME,
                               NUCLEOMORPH_GENE,
                               _5_METHYLCYTOSINE,
                               SL3_ACCEPTOR_SITE,
                               H4K4_TRIMETHYLATION_SITE,
                               INTRONIC_SPLICING_SILENCER,
                               N6_2_PRIME_O_DIMETHYLADENOSINE,
                               METAL_BINDING_SITE,
                               DMV5_MOTIF,
                               POLYPEPTIDE_BINDING_MOTIF,
                               GTT_TRINUCLEOTIDE_REPEAT,
                               U5_LTR_REGION,
                               N6_METHYLADENINE,
                               PYRIMIDINE_TRANSITION,
                               DH_REPEAT,
                               GLUTAMIC_ACID_TRNA_PRIMARY_TRANSCRIPT,
                               CONTIG_READ,
                               PRE_EDITED_REGION,
                               TRANSPOSABLE_ELEMENT_CDS,
                               H4R3_DIMETHYLATION_SITE,
                               CATMAT_LEFT_HANDED_FOUR,
                               ISOWYOSINE,
                               SS_RNA_VIRAL_SEQUENCE,
                               N6_N6_DIMETHYLADENOSINE,
                               VACUOLAR_SORTING_SIGNAL,
                               PROTEIN_CODING_PRIMARY_TRANSCRIPT,
                               SNRNA_PRIMARY_TRANSCRIPT,
                               TRANSCRIPTION_REGULATORY_REGION,
                               MIRNA_LOOP,
                               DOUBLE_STRANDED_DNA_CHROMOSOME,
                               INR1_MOTIF,
                               L_BOX,
                               N6_THREONYLCARBAMOYLADENOSINE,
                               SEVEN_AMINOMETHYL_SEVEN_DEAZAGUANOSINE,
                               ONE_METHYLGUANOSINE,
                               FIVE_PRIME_UTR,
                               UST_MATCH,
                               DCE_SII,
                               RR_TRACT,
                               EPIGENETICALLY_MODIFIED_REGION,
                               TELOMERASE_RNA_GENE,
                               PATERNALLY_IMPRINTED_GENE,
                               PSEUDOGENIC_TRNA,
                               PROPHAGE,
                               RNA_6S,
                               V_VDJ_C_CLUSTER,
                               INTERBAND,
                               INSERTION,
                               FIVE_ISOPENTENYLAMINOMETHYL_TWO_PRIME_O_METHYLURIDINE,
                               THREE_METHYLURIDINE,
                               H3K27_ACETYLATION_SITE,
                               CIRCULAR_DOUBLE_STRANDED_DNA_CHROMOSOME,
                               LINCRNA_GENE,
                               CHROMOSOME_BREAKAGE_SEQUENCE,
                               CTG_START_CODON,
                               H3K18_ACETYLATION_SITE,
                               POLYPEPTIDE,
                               H4K8_ACYLATION_SITE,
                               STOP_CODON_SIGNAL,
                               DELETION,
                               THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE,
                               SIGNATURE,
                               LEUCINE_TRNA_PRIMARY_TRANSCRIPT,
                               MT_GENE,
                               RESCUE_REGION,
                               NESTED_REPEAT,
                               MONOCISTRONIC_MRNA,
                               H3K36_DIMETHYLATION_SITE,
                               GENE,
                               MULTIPLEXING_SEQUENCE_IDENTIFIER,
                               MODIFIED_L_HISTIDINE,
                               X_ELEMENT,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA_ECF,
                               TRANSLATIONALLY_REGULATED_GENE,
                               METHYLWYOSINE,
                               ENGINEERED_PLASMID,
                               ENGINEERED_EPISOME,
                               SINGLE_STRAND_RESTRICTION_ENZYME_CLEAVAGE_SITE,
                               CDNA_CLONE,
                               MODIFIED_L_PHENYLALANINE,
                               CHROMOSOME_ARM,
                               TWO_METHYLTHIO_N6_THREONYL_CARBAMOYLADENOSINE,
                               NON_TRANSCRIBED_REGION,
                               UNEDITED_REGION,
                               BLUNT_END_RESTRICTION_ENZYME_CLEAVAGE_JUNCTION,
                               EXON_REGION,
                               RIBOSOME_ENTRY_SITE,
                               POLYPYRIMIDINE_TRACT,
                               SEQUENCE_LENGTH_VARIATION,
                               THREE_PRIME_RACE_CLONE,
                               GENOMIC_CLONE,
                               EDITED_CDS,
                               POLYA_SEQUENCE,
                               UNDERMODIFIED_HYDROXYWYBUTOSINE,
                               FRE,
                               ASX_TURN,
                               MIR_ENCODING_LNCRNA_PRIMARY_TRANSCRIPT,
                               THREE_PRIME_RST,
                               THREE_PRIME_D_SPACER,
                               CATALYTIC_RESIDUE,
                               SILENCER,
                               SL6_ACCEPTOR_SITE,
                               AMBISENSE_SSRNA_VIRAL_SEQUENCE,
                               TRANS_SPLICE_DONOR_SITE,
                               TMRNA_REGION,
                               BETA_TURN_TYPE_SIX_A_TWO,
                               ONE_TWO_PRIME_O_DIMETHYLGUANOSINE,
                               PROTEIN_STABILITY_ELEMENT,
                               NUCLEAR_LOCALIZATION_SIGNAL,
                               MATURE_TRANSCRIPT_REGION,
                               MRNA_WITH_PLUS_1_FRAMESHIFT,
                               ARS_CONSENSUS_SEQUENCE,
                               CIS_SPLICE_SITE,
                               EXPRESSED_SEQUENCE_ASSEMBLY,
                               EXON,
                               DISTAL_DUPLICATION,
                               MIRNA_STEM,
                               LNCRNA_PRIMARY_TRANSCRIPT,
                               ENZYMATIC_RNA,
                               TWO_PRIME_O_RIBOSYLADENOSINE_PHOSPHATE,
                               TRNA_GENE,
                               CAAT_SIGNAL,
                               POLYPEPTIDE_STRUCTURAL_REGION,
                               ATTENUATOR,
                               C_TERMINAL_REGION,
                               ORPHAN_CDS,
                               GRNA_GENE,
                               LINKAGE_GROUP,
                               SEQUENCE_UNCERTAINTY,
                               U14_SNORNA_PRIMARY_TRANSCRIPT,
                               SAGE_TAG,
                               TRANSLOCATION,
                               COMPLEX_SUBSTITUTION,
                               V_DJ_J_C_CLUSTER,
                               FIVE_ISOPENTENYLAMINOMETHYL_URIDINE,
                               LTR_COMPONENT,
                               DJ_GENE_SEGMENT,
                               THREONINE_TRNA_PRIMARY_TRANSCRIPT,
                               RST,
                               N6_METHYLADENOSINE,
                               N6_CIS_HYDROXYISOPENTENYL_ADENOSINE,
                               FIVE_HYDROXYMETHYLCYTIDINE,
                               NUCLEAR_GENE,
                               V_HEPTAMER,
                               BETA_TURN_RIGHT_HANDED_TYPE_TWO,
                               MODIFIED_DNA_BASE,
                               INR_MOTIF,
                               N4_ACETYLCYTIDINE,
                               PROPLASTID_GENE,
                               SEQUENCING_PRIMER,
                               ONE_METHYLPSEUDOURIDINE,
                               MODIFIED_AMINO_ACID_FEATURE,
                               BETA_BULGE_LOOP_FIVE,
                               HAMMERHEAD_RIBOZYME,
                               MCB,
                               NCRNA_GENE,
                               Z1_REGION,
                               H_PSEUDOKNOT,
                               OPEN_CHROMATIN_REGION,
                               MIR_ENCODING_VAULTRNA_PRIMARY_TRANSCRIPT,
                               TRANSPOSABLE_ELEMENT_PSEUDOGENE,
                               ZINC_REPRESSED_ELEMENT,
                               POSITIVELY_AUTOREGULATED_GENE,
                               D_J_C_CLUSTER,
                               LEUCOPLAST_GENE,
                               TMRNA_GENE,
                               INTRACHROMOSOMAL_BREAKPOINT,
                               ENGINEERED_TAG,
                               WYOSINE,
                               GENE_CASSETTE,
                               U2_SNRNA,
                               NON_LTR_RETROTRANSPOSON_POLYMERIC_TRACT,
                               POLYCISTRONIC_TRANSCRIPT,
                               BETA_TURN_TYPE_EIGHT,
                               MODIFIED_L_GLUTAMINE,
                               N6_GLYCINYLCARBAMOYLADENOSINE,
                               H2AK9_ACETYLATION_SITE,
                               LNA_OLIGO,
                               N_REGION,
                               CAGE_TAG,
                               THREE_PRIME_EST,
                               GENE_WITH_RECODED_MRNA,
                               BIOCHEMICAL_REGION_OF_PEPTIDE,
                               INTEGRATED_MOBILE_GENETIC_ELEMENT,
                               Y_RNA,
                               TRANS_SPLICED_TRANSCRIPT,
                               WC_BASE_PAIR,
                               ARCHAEOSINE,
                               PROMOTER,
                               ORIC,
                               GENE_REARRANGED_AT_DNA_LEVEL,
                               DDB_BOX,
                               H2BK120_ACETYLATION_SITE,
                               LEUCYL_TRNA,
                               TRNA_INTRON,
                               OPERON,
                               H3K4AC_ACETYLATION_SITE,
                               ER_RETENTION_SIGNAL,
                               FIVE_TWO_PRIME_O_DIMETHYLCYTIDINE,
                               TMRNA_ACCEPTOR_PIECE,
                               PIRNA_GENE,
                               POLYPEPTIDE_COPPER_ION_CONTACT_SITE,
                               GLUTAMYL_TRNA,
                               RNA_STABILITY_ELEMENT,
                               NUCLEASE_SENSITIVE_SITE,
                               CLIP,
                               _4_METHYLCYTOSINE,
                               ARGININE_TRNA_PRIMARY_TRANSCRIPT,
                               FLANKING_THREE_PRIME_QUADRUPLET_RECODING_SIGNAL,
                               MATURE_TRANSCRIPT,
                               TRANSLATED_NUCLEOTIDE_MATCH,
                               SUBSTITUTION,
                               GENE_FRAGMENT,
                               GENE_WITH_STOP_CODON_REDEFINED_AS_PYRROLYSINE,
                               MITOCHONDRIAL_CHROMOSOME,
                               MODIFIED_L_ALANINE,
                               REGIONAL_CENTROMERE_CENTRAL_CORE,
                               MATING_TYPE_REGION,
                               NUCLEOTIDE_TO_PROTEIN_BINDING_SITE,
                               MODIFIED_URIDINE,
                               RFLP_FRAGMENT,
                               CYTOPLASMIC_POLYPEPTIDE_REGION,
                               D_DJ_C_CLUSTER,
                               TRANSMEMBRANE_HELIX,
                               V_REGION,
                               GAP,
                               V_J_CLUSTER,
                               PCB,
                               INTEGRATION_EXCISION_SITE,
                               LOCUS_CONTROL_REGION,
                               C_REGION,
                               DNA_CONSTRAINT_SEQUENCE,
                               LONG_TERMINAL_REPEAT,
                               MATRIX_ATTACHMENT_SITE,
                               TWO_THIOCYTIDINE,
                               REVERSE_PRIMER,
                               GLYCINE,
                               TELOMERIC_REPEAT,
                               THREE_PRIME_CODING_EXON_CODING_REGION,
                               MINICIRCLE_GENE,
                               POLYPEPTIDE_NEST_MOTIF,
                               NON_CANONICAL_THREE_PRIME_SPLICE_SITE,
                               CROSS_GENOME_MATCH,
                               ACE2_UAS,
                               CDRE_MOTIF,
                               FIVE_ISOPENTENYLAMINOMETHYL_TWO_THIOURIDINE,
                               POST_TRANSLATIONALLY_REGULATED_GENE,
                               PLASMID_VECTOR,
                               VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_SPACER,
                               TRANSCRIBED_FRAGMENT,
                               GENE_WITH_EDITED_TRANSCRIPT,
                               RECOMBINATIONALLY_REARRANGED_GENE,
                               ORIGIN_OF_REPLICATION,
                               DMV1_MOTIF,
                               MINUS_10_SIGNAL,
                               VJ_GENE_SEGMENT,
                               PROCESSED_TRANSCRIPT,
                               NUCLEASE_HYPERSENSITIVE_SITE,
                               TCT_MOTIF,
                               SEVEN_CYANO_SEVEN_DEAZAGUANOSINE,
                               TF_BINDING_SITE,
                               MNV,
                               V_VJ_J_C_CLUSTER,
                               RECOMBINATION_FEATURE_OF_REARRANGED_GENE,
                               GATA_BOX,
                               RNAI_REAGENT,
                               REPEAT_UNIT,
                               SORF,
                               TNA_OLIGO,
                               PLASMID,
                               H3K27_MONOMETHYLATION_SITE,
                               REGIONAL_CENTROMERE_OUTER_REPEAT_REGION,
                               MAJOR_TSS,
                               MODIFIED_ADENOSINE,
                               GENE_WITH_STOP_CODON_REDEFINED_AS_SELENOCYSTEINE,
                               FORWARD_PRIMER,
                               UNASSIGNED_SUPERCONTIG,
                               GUIDE_RNA_REGION,
                               READ,
                               VDJ_C_CLUSTER,
                               IMMATURE_PEPTIDE_REGION,
                               TATA_BOX,
                               BAC,
                               ASPARAGINE,
                               LIGATION_BASED_READ,
                               NONAMER_OF_RECOMBINATION_FEATURE_OF_VERTEBRATE_IMMUNE_SYSTEM_GENE,
                               H4K20_MONOMETHYLATION_SITE,
                               PAC_END,
                               MANNOSYL_QUEUOSINE,
                               UORF,
                               _5_CARBOXYLCYTOSINE,
                               POLYPEPTIDE_METAL_CONTACT,
                               CDS_FRAGMENT,
                               BETA_BULGE_LOOP_SIX,
                               MIRNA_TARGET_SITE,
                               TRANS_SPLICE_SITE,
                               RNA_MOTIF,
                               MODIFIED_L_THREONINE,
                               V_D_J_CLUSTER,
                               SATELLITE_DNA,
                               MODIFIED_INOSINE,
                               PSEUDOURIDINE,
                               VDJ_J_CLUSTER,
                               GAGA_MOTIF,
                               PRIMER_BINDING_SITE,
                               N6_METHYL_N6_THREONYLCARBAMOYLADENOSINE,
                               RRNA_28S,
                               DNA_MOTIF,
                               MICROARRAY_OLIGO,
                               TRANSVERSION,
                               DNA_CHROMOSOME,
                               INSULATOR_BINDING_SITE,
                               RECODING_STIMULATORY_REGION,
                               POSSIBLE_ASSEMBLY_ERROR,
                               DISTANT_THREE_PRIME_RECODING_SIGNAL,
                               PROVIRAL_GENE,
                               PHAGEMID,
                               MACRONUCLEAR_CHROMOSOME,
                               V_D_DJ_J_C_CLUSTER,
                               MOR,
                               TMRNA_PRIMARY_TRANSCRIPT,
                               THREE_PRIME_CIS_SPLICE_SITE,
                               U6_SNRNA,
                               SEQUENCE_FEATURE,
                               DNAZYME,
                               ALPHA_HELIX,
                               HYDROPHOBIC_REGION_OF_PEPTIDE,
                               RNA_HOOK_TURN,
                               LNC_RNA,
                               STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE,
                               GENE_MEMBER_REGION,
                               SERINE_TRNA_PRIMARY_TRANSCRIPT,
                               CORE_PROMOTER_ELEMENT,
                               PSEUDOGENIC_RRNA,
                               RNA_POLYMERASE_III_TATA_BOX,
                               DJ_J_CLUSTER,
                               POLYPEPTIDE_TUNGSTEN_ION_CONTACT_SITE,
                               SS_OLIGO,
                               U11_SNRNA,
                               VJ_C_CLUSTER,
                               ACTIVE_PEPTIDE,
                               CYSTEINYL_TRNA,
                               ST_TURN_RIGHT_HANDED_TYPE_ONE,
                               CONFORMATIONAL_SWITCH,
                               SELENOCYSTEINE_TRNA_PRIMARY_TRANSCRIPT,
                               PEPTIDE_LOCALIZATION_SIGNAL,
                               THREE_PRIME_LTR_COMPONENT,
                               BIOMATERIAL_REGION,
                               GENE_WITH_NON_CANONICAL_START_CODON,
                               READ_PAIR,
                               DINUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE,
                               INVERTED_TANDEM_DUPLICATION,
                               GOLDEN_PATH,
                               TRANSIT_PEPTIDE,
                               NOVEL_SEQUENCE_INSERTION,
                               RRNA_PRIMARY_TRANSCRIPT,
                               CDS,
                               MATURE_PROTEIN_REGION,
                               EPIGENETICALLY_MODIFIED_GENE,
                               H3K9_METHYLATION_SITE,
                               FOUR_BP_START_CODON,
                               GENE_WITH_START_CODON_CUG,
                               MUTATED_VARIANT_SITE,
                               MINUS_35_SIGNAL,
                               ARS,
                               ASSEMBLY_COMPONENT,
                               TRANSPOSABLE_ELEMENT,
                               J_HEPTAMER,
                               CT_GENE,
                               NUCLEAR_CHROMOSOME,
                               H2BK15_ACETYLATION_SITE,
                               SERINE_THREONINE_STAPLE_MOTIF,
                               STS_MAP,
                               CYANELLE_CHROMOSOME,
                               TWINTRON,
                               CRYPTOGENE,
                               KINETOPLAST_GENE,
                               DS_DNA_VIRAL_SEQUENCE,
                               HISTONE_BINDING_SITE,
                               BACTERIAL_TERMINATOR,
                               POLYMERASE_SYNTHESIS_READ,
                               TWO_METHYLTHIO_N6_METHYLADENOSINE,
                               U2_INTRON,
                               FUNCTIONAL_CANDIDATE_GENE,
                               HISTONE_UBIQITINATION_SITE,
                               THREE_PRIME_D_NONAMER,
                               NCRNA,
                               H2BK12_ACETYLATION_SITE,
                               FIVE_CARBOXYHYDROXYMETHYL_URIDINE,
                               RETROTRANSPOSON,
                               INTRONIC_SPLICE_ENHANCER,
                               SPECIFIC_RECOMBINATION_SITE,
                               CANDIDATE_GENE,
                               RRE_RNA,
                               M26_BINDING_SITE,
                               ISOLEUCYL_TRNA,
                               A_TO_T_TRANSVERSION,
                               HISTIDINE_TRNA_PRIMARY_TRANSCRIPT,
                               MINISATELLITE,
                               FIVE_CARBOXYMETHYLAMINOMETHYL_TWO_PRIME_O_METHYLURIDINE,
                               SPOT_42_RNA,
                               H3K79_DIMETHYLATION_SITE,
                               DISTAL_PROMOTER_ELEMENT,
                               MODIFIED_L_SERINE,
                               POLYPEPTIDE_CALCIUM_ION_CONTACT_SITE,
                               CCAAT_MOTIF,
                               METABOLIC_ISLAND,
                               VAULTRNA_PRIMARY_TRANSCRIPT,
                               RNAPOL_III_PROMOTER_TYPE_1,
                               PSEUDOGENIC_GENE_SEGMENT,
                               WILD_TYPE_RESCUE_GENE,
                               GROUP_1_INTRON_HOMING_ENDONUCLEASE_TARGET_REGION,
                               RNASE_P_RNA_GENE,
                               V_DJ_CLUSTER,
                               RNASE_P_RNA,
                               CLEAVED_FOR_GPI_ANCHOR_REGION,
                               ORF,
                               RESCUE_MINI_GENE,
                               HERITABLE_PHENOTYPIC_MARKER,
                               BAC_CLONED_GENOMIC_INSERT,
                               POLYPEPTIDE_NICKEL_ION_CONTACT_SITE,
                               HAPLOTYPE_BLOCK,
                               THREE_PRIME_INTRON,
                               INTERMEDIATE_ELEMENT,
                               DESTRUCTION_BOX,
                               H3K27_DIMETHYLATION_SITE,
                               S_GNA_OLIGO,
                               GENE_COMPONENT_REGION,
                               V_VJ_CLUSTER,
                               ORIV,
                               SNRNA_GENE,
                               FOUR_DEMETHYLWYOSINE,
                               TELOMERE,
                               EXTRAMEMBRANE_POLYPEPTIDE_REGION,
                               RESTRICTION_ENZYME_CLEAVAGE_JUNCTION,
                               BIDIRECTIONAL_PROMOTER,
                               UPD,
                               SYNTENIC_REGION,
                               ANTIPARALLEL_BETA_STRAND,
                               ENHANCER_BOUND_BY_FACTOR,
                               DPE1_MOTIF,
                               FLOXED_GENE,
                               ALANINE_TRNA_PRIMARY_TRANSCRIPT,
                               NESTED_TRANSPOSON,
                               SHRNA,
                               RIBOTHYMIDINE,
                               CLONED_GENOMIC_INSERT,
                               EXPERIMENTAL_FEATURE,
                               PROTEIN_BINDING_SITE,
                               SCHELLMANN_LOOP_SIX,
                               ENDOSOMAL_LOCALIZATION_SIGNAL,
                               _8_OXOGUANINE,
                               GLYCYL_TRNA,
                               INTROGRESSED_CHROMOSOME_REGION,
                               ONE_METHYLADENOSINE,
                               FIVE_METHYLAMINOMETHYLURIDINE,
                               BETA_STRAND,
                               HOMING_ENDONUCLEASE_BINDING_SITE,
                               TRYPTOPHANYL_TRNA,
                               CYANELLE_GENE,
                               MODIFIED_GUANINE,
                               GALACTOSYL_QUEUOSINE,
                               MATCH,
                               POINT_MUTATION,
                               HISTONE_4_ACYLATION_SITE,
                               CDS_PREDICTED,
                               TRANSCRIPT_BOUND_BY_NUCLEIC_ACID,
                               BETA_TURN_TYPE_SIX_A_ONE,
                               TRANSGENIC_TRANSPOSABLE_ELEMENT,
                               PYROSEQUENCED_READ,
                               ONE_METHYL_THREE_THREE_AMINO_THREE_CARBOXYPROPYL_PSEUDOURIDINE,
                               NUCLEOTIDE_MOTIF,
                               V_VDJ_J_C_CLUSTER,
                               ALANINE,
                               H4K5_ACYLATION_SITE,
                               RECURSIVE_SPLICE_SITE,
                               FOREIGN_GENE,
                               VJ_J_C_CLUSTER,
                               FIVE_METHOXYCARBONYLMETHYL_TWO_THIOURIDINE,
                               H3K9_TRIMETHYLATION_SITE,
                               POLYPEPTIDE_STRUCTURAL_MOTIF,
                               MATING_TYPE_REGION_REPLICATION_FORK_BARRIER,
                               FIVE_FORMYL_TWO_PRIME_O_METHYLCYTIDINE,
                               COPY_NUMBER_GAIN,
                               BETA_TURN_TYPE_SIX,
                               RRNA_5_8S,
                               TWO_PRIME_O_METHYLGUANOSINE,
                               DICISTRONIC_TRANSCRIPT,
                               THREE_PRIME_CODING_EXON_NONCODING_REGION,
                               GROUP_IIA_INTRON,
                               CAPPED_PRIMARY_TRANSCRIPT,
                               FIVE_PRIME_CODING_EXON_NONCODING_REGION,
                               RRNA_CLEAVAGE_SNORNA_PRIMARY_TRANSCRIPT,
                               A_BOX_TYPE_2,
                               PSE_MOTIF,
                               ENGINEERED_GENE,
                               REPEAT_REGION,
                               BAC_END,
                               LINE_ELEMENT,
                               FOREIGN_TRANSPOSABLE_ELEMENT,
                               OCTAMER_MOTIF,
                               PURINE_TO_PYRIMIDINE_TRANSVERSION,
                               VJ_J_CLUSTER,
                               CENTROMERIC_REPEAT,
                               MEMBRANE_PEPTIDE_LOOP,
                               SRP_RNA_PRIMARY_TRANSCRIPT,
                               H3K36_TRIMETHYLATION_SITE,
                               CHROMOSOME_PART,
                               GCVB_RNA,
                               BETA_TURN_TYPE_SIX_A,
                               THREE_METHYLPSEUDOURIDINE,
                               TILING_PATH,
                               PRIMARY_TRANSCRIPT,
                               PROCESSED_PSEUDOGENE,
                               FIVE_PRIME_CODING_EXON,
                               GOLDEN_PATH_FRAGMENT,
                               EXPERIMENTAL_RESULT_REGION,
                               CATMAT_LEFT_HANDED_THREE,
                               GENOMIC_ISLAND,
                               TRANSCRIPTIONAL_CIS_REGULATORY_REGION,
                               QUEUOSINE,
                               BASE_PAIR,
                               PRE_MIRNA,
                               SERINE_THREONINE_MOTIF,
                               RRNA_GENE,
                               SHINE_DALGARNO_SEQUENCE,
                               SILENCED_GENE,
                               VDJ_GENE_SEGMENT,
                               INVERSION_SITE_PART,
                               RNA_POLYMERASE_II_TATA_BOX,
                               SGRNA,
                               INTERIOR_INTRON,
                               SNV,
                               PREDICTED_GENE,
                               N6_ISOPENTENYLADENOSINE,
                               TERRA,
                               INTRONIC_REGULATORY_REGION,
                               ASYMMETRIC_RNA_INTERNAL_LOOP,
                               SL11_ACCEPTOR_SITE,
                               DSRA_RNA,
                               DIF_SITE,
                               RRNA_18S,
                               MIR_ENCODING_SNORNA_PRIMARY_TRANSCRIPT,
                               THREE_PRIME_REPEAT_RECODING_SIGNAL,
                               EDITED_MRNA,
                               MOBILE_INTRON,
                               TERMINAL_INVERTED_REPEAT_ELEMENT,
                               ISOLEUCINE,
                               TYROSINE,
                               MINI_GENE,
                               DCAPS_PRIMER,
                               DS_RNA_VIRAL_SEQUENCE,
                               H3K23_DIMETHYLATION_SITE,
                               NO_OUTPUT,
                               TRANSCRIBED_CLUSTER,
                               TWO_METHYLADENOSINE,
                               N2_2_PRIME_O_DIMETHYLGUANOSINE,
                               RESTRICTION_ENZYME_ASSEMBLY_SCAR,
                               RETINOIC_ACID_RESPONSIVE_ELEMENT,
                               GENE_SILENCED_BY_DNA_METHYLATION,
                               OUTRON,
                               SNP,
                               TWO_THIO_TWO_PRIME_O_METHYLURIDINE,
                               SEQUENCE_REARRANGEMENT_FEATURE,
                               DIHYDROURIDINE,
                               KOZAK_SEQUENCE,
                               MOLECULAR_CONTACT_REGION,
                               H3K79_TRIMETHYLATION_SITE,
                               FIVE_METHOXYCARBONYLMETHYL_TWO_PRIME_O_METHYLURIDINE,
                               GENE_CASSETTE_ARRAY,
                               CAGE_CLUSTER,
                               FIVE_PRIME_RECODING_SITE,
                               VIRAL_SEQUENCE,
                               GLUTAMIC_ACID,
                               EXON_OF_SINGLE_EXON_GENE,
                               D_LOOP,
                               MIRNA_PRIMARY_TRANSCRIPT_REGION,
                               N_TERMINAL_REGION,
                               MITOCHONDRIAL_SUPERCONTIG,
                               N4_2_PRIME_O_DIMETHYLCYTIDINE,
                               CLASS_I_RNA,
                               CLEAVED_PEPTIDE_REGION,
                               FIVE_PRIME_RESTRICTION_ENZYME_JUNCTION,
                               LINEAR_DOUBLE_STRANDED_DNA_CHROMOSOME,
                               V_VJ_J_CLUSTER,
                               TETRALOOP,
                               MOBILE_ELEMENT_INSERTION,
                               REPLICATION_REGULATORY_REGION,
                               D_GENE_RECOMBINATION_FEATURE,
                               TWO_PRIME_O_RIBOSYLGUANOSINE_PHOSPHATE,
                               C_CLUSTER,
                               MICRONUCLEAR_CHROMOSOME,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA54,
                               PEPTIDE_COIL,
                               NC_CONSERVED_REGION,
                               RNA_SEQUENCE_SECONDARY_STRUCTURE,
                               RNASE_MRP_RNA_GENE,
                               REGION,
                               PLUS_2_TRANSLATIONAL_FRAMESHIFT,
                               J_CLUSTER,
                               SEQUENCE_ASSEMBLY,
                               MRNA,
                               NONCODING_REGION_OF_EXON,
                               DRE_MOTIF,
                               CENTRAL_HYDROPHOBIC_REGION_OF_SIGNAL_PEPTIDE,
                               BRUNO_RESPONSE_ELEMENT,
                               START_CODON,
                               HSE,
                               D_J_CLUSTER,
                               TRYPTOPHAN_TRNA_PRIMARY_TRANSCRIPT,
                               HOMOLOGOUS_REGION,
                               FRT_SITE,
                               WOBBLE_BASE_PAIR,
                               CLONED_CDNA_INSERT,
                               MODIFIED_GUANOSINE,
                               CLONE_INSERT,
                               QTL,
                               ARRET,
                               PROPEPTIDE,
                               CODING_CONSERVED_REGION,
                               RRNA,
                               VERTEBRATE_IMMUNE_SYSTEM_GENE_RECOMBINATION_FEATURE,
                               CYSTEINE_TRNA_PRIMARY_TRANSCRIPT,
                               FIVE_TAURINOMETHYLURIDINE,
                               G_TO_C_TRANSVERSION,
                               RESCUE_GENE,
                               NDM2_MOTIF,
                               SONICATE_FRAGMENT,
                               ENGINEERED_RESCUE_REGION,
                               INSERTION_BREAKPOINT,
                               GT_DINUCLEOTIDE_REPEAT,
                               C_BOX,
                               TWO_PRIME_O_METHYLADENOSINE,
                               H3K56_ACETYLATION_SITE,
                               MODIFIED_L_SELENOCYSTEINE,
                               FIVE_FORMYLCYTIDINE,
                               INSULATOR,
                               BETA_TURN_TYPE_SIX_B,
                               POLYPEPTIDE_CATALYTIC_MOTIF,
                               H3K27_ACYLATION_SITE,
                               H3R2_DIMETHYLATION_SITE,
                               NUCLEOTIDE_MATCH,
                               R_LTR_REGION,
                               SEQUENCE_SECONDARY_STRUCTURE,
                               INVERSION_SITE,
                               SPLICE_REGION,
                               SCRNA_PRIMARY_TRANSCRIPT,
                               TERMINATOR,
                               POLINTON,
                               INTRON_DOMAIN,
                               FIVE_PRIME_NONCODING_EXON,
                               TASIRNA,
                               LYSYL_TRNA,
                               LEUCOPLAST_CHROMOSOME,
                               GENE_WITH_TRANS_SPLICED_TRANSCRIPT,
                               TANDEM_REPEAT,
                               N2_N2_2_PRIME_O_TRIMETHYLGUANOSINE,
                               ATTC_SITE,
                               RNAPOL_II_PROMOTER,
                               MATHEMATICALLY_DEFINED_REPEAT,
                               DNA_SEQUENCE_SECONDARY_STRUCTURE,
                               SPLICEOSOMAL_INTRON,
                               A_TO_G_TRANSITION,
                               H2BK5_MONOMETHYLATION_SITE,
                               MIR_ENCODING_Y_RNA_PRIMARY_TRANSCRIPT,
                               V_D_J_C_CLUSTER,
                               GUIDE_RNA,
                               EDITED_TRANSCRIPT_BY_A_TO_I_SUBSTITUTION,
                               IRON_RESPONSIVE_ELEMENT,
                               V_SPACER,
                               NUCLEAR_MT_PSEUDOGENE,
                               VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_REARRANGED_GENE_CLUSTER,
                               CHLOROPLAST_DNA_READ,
                               TMRNA_CODING_PIECE,
                               BLOCKED_READING_FRAME,
                               SILENT_MATING_TYPE_CASSETTE_ARRAY,
                               SCHELLMANN_LOOP_SEVEN,
                               STOP_CODON,
                               POLYPEPTIDE_NEST_RIGHT_LEFT_MOTIF,
                               BRED_MOTIF,
                               THREE_PRIME_CLIP,
                               ORIT,
                               MRNA_CONTIG,
                               LARGE_SUBUNIT_RRNA,
                               NUCLEAR_EXPORT_SIGNAL,
                               INVERTED_REPEAT,
                               LTR_RETROTRANSPOSON,
                               TELOMERASE_RNA,
                               G_TO_T_TRANSVERSION,
                               RH_MAP,
                               T_LOOP,
                               MRNA_WITH_FRAMESHIFT,
                               SOLO_LTR,
                               BOUNDARY_ELEMENT,
                               NUCLEAR_RIM_LOCALIZATION_SIGNAL,
                               ENGINEERED_TRANSPOSABLE_ELEMENT,
                               ARGINYL_TRNA,
                               TRANSCRIPT,
                               MATING_TYPE_REGION_MOTIF,
                               MODIFIED_L_CYSTEINE,
                               MONOMERIC_REPEAT,
                               INTERNAL_TRANSCRIBED_SPACER_REGION,
                               METHYLINOSINE,
                               FLANKING_REGION,
                               PI_HELIX,
                               ASPARAGINE_TRNA_PRIMARY_TRANSCRIPT,
                               BINDING_SITE,
                               X_ELEMENT_COMBINATORIAL_REPEAT,
                               ALLELICALLY_EXCLUDED_GENE,
                               MATERNALLY_IMPRINTED_GENE,
                               FIVE_AMINOMETHYL_TWO_THIOURIDINE,
                               MOBILE_GENETIC_ELEMENT,
                               AACCCT_BOX,
                               VALYL_TRNA,
                               EXEMPLAR_MRNA,
                               ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT_GENE,
                               TRANSITION,
                               DISPERSED_REPEAT,
                               CASSETTE_PSEUDOGENE,
                               EXONIC_SPLICE_ENHANCER,
                               PYRROLYSINE,
                               FIVE_METHOXYURIDINE,
                               MRNA_WITH_MINUS_1_FRAMESHIFT,
                               H3K9_ACETYLATION_SITE,
                               ANDROGEN_RESPONSE_ELEMENT,
                               APTAMER,
                               FIVE_PRIME_D_NONAMER,
                               GENE_SILENCED_BY_HISTONE_METHYLATION,
                               OPERATOR,
                               DCE,
                               CONSENSUS_GDNA,
                               CDS_SUPPORTED_BY_EST_OR_CDNA_DATA,
                               FIVE_TAURINOMETHYL_TWO_THIOURIDINE,
                               SNORNA_PRIMARY_TRANSCRIPT,
                               RDNA_REPLICATION_FORK_BARRIER,
                               VDJ_J_C_CLUSTER,
                               H3R2_MONOMETHYLATION_SITE,
                               STOP_CODON_REDEFINED_AS_PYRROLYSINE,
                               EDITED_TRANSCRIPT_FEATURE,
                               SL4_ACCEPTOR_SITE,
                               THREE_THREE_AMINO_THREE_CARBOXYPROPYL_URIDINE,
                               CATMAT_RIGHT_HANDED_FOUR,
                               U5_FIVE_PRIME_LTR_REGION,
                               TRNA_PRIMARY_TRANSCRIPT,
                               DNA_BINDING_SITE,
                               FINGERPRINT_MAP,
                               GENE_WITH_DICISTRONIC_TRANSCRIPT,
                               UGA_STOP_CODON_SIGNAL,
                               H3K23_ACYLATION_SITE,
                               THREE_PRIME_RECODING_SITE,
                               INTERIOR_EXON,
                               THREE_PRIME_UTR_INTRON,
                               MAXICIRCLE_GENE,
                               N2_N2_DIMETHYLGUANOSINE,
                               T3_RNA_POLYMERASE_PROMOTER,
                               J_SPACER,
                               SITE_SPECIFIC_RECOMBINATION_TARGET_REGION,
                               RECODING_PSEUDOKNOT,
                               GENOMIC_DNA_READ,
                               MICF_RNA,
                               IMMUNOGLOBULIN_REGION,
                               SPLICEOSOMAL_INTRON_REGION,
                               MRNA_RECODED_BY_CODON_REDEFINITION,
                               RASIRNA,
                               TANDEM_DUPLICATION,
                               POLYPEPTIDE_COBALT_ION_CONTACT_SITE,
                               INTRINSICALLY_UNSTRUCTURED_POLYPEPTIDE_REGION,
                               CENTROMERE_DNA_ELEMENT_II,
                               SEQUENCE_ALTERATION,
                               VECTOR_REPLICON,
                               FIVE_PRIME_EST,
                               AMPLIFICATION_ORIGIN,
                               CONSENSUS_REGION,
                               N6_N6_2_PRIME_O_TRIMETHYLADENOSINE,
                               A_TO_C_TRANSVERSION,
                               CHROMOPLAST_GENE,
                               COILED_COIL,
                               TRYPTOPHAN,
                               V_D_DJ_C_CLUSTER,
                               ASX_TURN_LEFT_HANDED_TYPE_TWO,
                               X_REGION,
                               TRANSPOSON_FRAGMENT,
                               NONCODING_EXON,
                               SL2_ACCEPTOR_SITE,
                               HISTONE_MODIFICATION,
                               U_BOX,
                               DOUBLE_STRANDED_RNA_CHROMOSOME,
                               PYRROLYSYL_TRNA,
                               NUCLEASE_BINDING_SITE,
                               MIRNA,
                               FIVE_PRIME_FLANKING_REGION,
                               CRISPR,
                               SARCIN_LIKE_RNA_MOTIF,
                               RHO_DEPENDENT_BACTERIAL_TERMINATOR,
                               H3K4_MONOMETHYLATION_SITE,
                               RNA_POLYMERASE_PROMOTER,
                               SINGLE_STRANDED_RNA_CHROMOSOME,
                               N4_METHYLCYTIDINE,
                               THREE_PRIME_RESTRICTION_ENZYME_JUNCTION,
                               POLYPEPTIDE_TURN_MOTIF,
                               THREE_PRIME_UST,
                               DCE_SI,
                               ST_TURN_LEFT_HANDED_TYPE_ONE,
                               HIGH_IDENTITY_REGION,
                               H3K20_TRIMETHYLATION_SITE,
                               MODIFIED_ADENINE,
                               SEQUENCE_CONFLICT,
                               AMINO_ACID,
                               ANTICODON,
                               SIMPLE_SEQUENCE_LENGTH_VARIATION,
                               BAC_READ_CONTIG,
                               SIGNAL_ANCHOR,
                               INTERNAL_RIBOSOME_ENTRY_SITE,
                               V_D_DJ_J_CLUSTER,
                               SRP_RNA_GENE,
                               VERTEBRATE_IMMUNOGLOBULIN_T_CELL_RECEPTOR_SEGMENT,
                               INDUCIBLE_PROMOTER,
                               V_DJ_J_CLUSTER,
                               DUPLICATED_PSEUDOGENE,
                               _5_HYDROXYMETHYLCYTOSINE,
                               REGIONAL_CENTROMERE_INNER_REPEAT_REGION,
                               FIVE_PRIME_LTR,
                               U4ATAC_SNRNA,
                               MATING_TYPE_M_BOX,
                               EXONIC_SPLICING_SILENCER,
                               ANCHOR_REGION,
                               CONSERVED_REGION,
                               MAXICIRCLE,
                               TNAORF,
                               GENE_WITH_POLYADENYLATED_MRNA,
                               SMALL_SUBUNIT_RRNA,
                               FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE,
                               SL8_ACCEPTOR_SITE,
                               REMARK,
                               SCRNA_GENE,
                               V_CLUSTER,
                               INTRON,
                               INTEIN_ENCODING_REGION,
                               OVERLAPPING_EST_SET,
                               CURE,
                               GENE_SILENCED_BY_HISTONE_MODIFICATION,
                               URIDINE_FIVE_OXYACETIC_ACID_METHYL_ESTER,
                               LYSINE_TRNA_PRIMARY_TRANSCRIPT,
                               TRANS_SPLICED_MRNA,
                               RHO_INDEPENDENT_BACTERIAL_TERMINATOR,
                               MODIFIED_L_ARGININE,
                               GENE_ARRAY,
                               MODIFIED_L_ISOLEUCINE,
                               T_TO_C_TRANSITION,
                               N4_ACETYL_2_PRIME_O_METHYLCYTIDINE,
                               TWO_PRIME_O_METHYLCYTIDINE,
                               TRNA,
                               VALINE,
                               DELETION_JUNCTION,
                               PSEUDOGENIC_TRANSCRIPT,
                               COPY_NUMBER_LOSS,
                               U4_SNRNA,
                               VALINE_TRNA_PRIMARY_TRANSCRIPT,
                               COSMID,
                               THREE_TEN_HELIX,
                               EUKARYOTIC_TERMINATOR,
                               FIVE_METHYL_2_THIOURIDINE,
                               FIVE_CARBOXYHYDROXYMETHYL_URIDINE_METHYL_ESTER,
                               CONJUGATIVE_TRANSPOSON,
                               PRIMER,
                               EPISOME,
                               ULTRACONTIG,
                               DJ_J_C_CLUSTER,
                               COMPLEX_STRUCTURAL_ALTERATION,
                               POLYPEPTIDE_ZINC_ION_CONTACT_SITE,
                               CODING_END,
                               BRANCH_SITE,
                               TWO_PRIME_O_METHYLINOSINE,
                               METHIONINE,
                               UTR_INTRON,
                               PROMOTER_ELEMENT,
                               POST_TRANSLATIONALLY_MODIFIED_REGION,
                               D_DJ_J_CLUSTER,
                               FOLDBACK_ELEMENT,
                               G_TO_A_TRANSITION,
                               GLUTAMINYL_TRNA,
                               FIVE_CARBOXYMETHYLURIDINE,
                               GENE_WITH_DICISTRONIC_MRNA,
                               ATTP_SITE,
                               TRANSLOCATION_BREAKPOINT,
                               REGULATORY_PROMOTER_ELEMENT,
                               MATERNAL_UNIPARENTAL_DISOMY,
                               K_TURN_RNA_MOTIF,
                               NON_CYTOPLASMIC_POLYPEPTIDE_REGION,
                               CIRCULAR_DOUBLE_STRANDED_RNA_CHROMOSOME,
                               GAMMA_TURN_CLASSIC,
                               RNAPOL_III_PROMOTER,
                               PSEUDOGENIC_EXON,
                               NATURAL_PLASMID,
                               APICOPLAST_GENE,
                               EDITING_DOMAIN,
                               INTRON_BASE_5,
                               MODIFIED_RNA_BASE_FEATURE,
                               LEUCINE,
                               TRANSGENE,
                               THREE_PRIME_FLANKING_REGION,
                               H3K4_METHYLATION_SITE,
                               INTERCHROMOSOMAL_BREAKPOINT,
                               PHOSPHORYLATION_SITE,
                               REGIONAL_CENTROMERE,
                               RAPD,
                               H3K14_ACETYLATION_SITE,
                               POLYPEPTIDE_DOMAIN,
                               HOOGSTEEN_BASE_PAIR,
                               PROXIMAL_PROMOTER_ELEMENT,
                               C_D_BOX_SNORNA,
                               RESTRICTION_ENZYME_RECOGNITION_SITE,
                               POLYPEPTIDE_MOTIF,
                               CSRB_RSMB_RNA,
                               TMRNA,
                               POSITIONAL_CANDIDATE_GENE,
                               LARIAT_INTRON,
                               LYSINE,
                               ATTR_SITE,
                               SUGAR_EDGE_BASE_PAIR,
                               BETA_TURN_LEFT_HANDED_TYPE_ONE,
                               CHIP_SEQ_REGION,
                               TRANSGENIC_INSERTION,
                               SIGNAL_PEPTIDE,
                               W_REGION,
                               S_REGION,
                               RRNA_21S,
                               FUSION_GENE,
                               CHROMOSOMAL_STRUCTURAL_ELEMENT,
                               N6_HYDROXYNORVALYLCARBAMOYLADENOSINE,
                               H2B_UBIQUITINATION_SITE,
                               RST_MATCH,
                               CCA_TAIL,
                               SCHELLMANN_LOOP,
                               GENE_SILENCED_BY_DNA_MODIFICATION,
                               ENGINEERED_INSERT,
                               UNITARY_PSEUDOGENE,
                               FIVE_PRIME_RST,
                               U12_SNRNA,
                               U14_SNORNA,
                               MODIFIED_L_LEUCINE,
                               RESOLUTION_SITE,
                               SL7_ACCEPTOR_SITE,
                               PROBE,
                               POLYPEPTIDE_SECONDARY_STRUCTURE,
                               POLYPEPTIDE_SEQUENCING_INFORMATION,
                               V_GENE_RECOMBINATION_FEATURE,
                               FIVE_CARBAMOYLMETHYL_TWO_PRIME_O_METHYLURIDINE,
                               TRINUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE,
                               CLONED_REGION,
                               MIR_ENCODING_SHRNA_PRIMARY_TRANSCRIPT,
                               GENE_SILENCED_BY_HISTONE_DEACETYLATION,
                               CIRCULAR_SINGLE_STRANDED_RNA_CHROMOSOME,
                               H3K9_MONOMETHYLATION_SITE,
                               ENDOGENOUS_RETROVIRAL_GENE,
                               FIVE_METHYLCYTIDINE,
                               STEM_LOOP,
                               PROLYL_TRNA,
                               SCRNA,
                               ORTHOLOGOUS_REGION,
                               FOSMID,
                               INTEIN,
                               GENE_WITH_MRNA_RECODED_BY_TRANSLATIONAL_BYPASS,
                               IRON_REPRESSED_GATA_ELEMENT,
                               HISTIDINE,
                               POLYPEPTIDE_MAGNESIUM_ION_CONTACT_SITE,
                               OKAZAKI_FRAGMENT,
                               UST,
                               H4K12_ACETYLATION_SITE,
                               POLYPEPTIDE_CONSERVED_MOTIF,
                               NATURAL_VARIANT_SITE,
                               H3K4_DIMETHYLATION_SITE,
                               OVERLAPPING_FEATURE_SET,
                               GLYCINE_TRNA_PRIMARY_TRANSCRIPT,
                               RIBONUCLEASE_SITE,
                               PLASTID_GENE,
                               VD_GENE_SEGMENT,
                               V_VDJ_J_CLUSTER,
                               SINGLE_STRANDED_DNA_CHROMOSOME,
                               POLYPEPTIDE_NEST_LEFT_RIGHT_MOTIF,
                               ASX_MOTIF,
                               RNAPOL_I_PROMOTER,
                               ASSEMBLY_ERROR_CORRECTION,
                               PHENYLALANYL_TRNA,
                               BLUNT_END_RESTRICTION_ENZYME_CLEAVAGE_SITE,
                               UAG_STOP_CODON_SIGNAL,
                               POLYPEPTIDE_REGION,
                               RESTRICTION_ENZYME_THREE_PRIME_SINGLE_STRAND_OVERHANG,
                               ST_TURN_RIGHT_HANDED_TYPE_TWO,
                               RNA_INTERNAL_LOOP,
                               POLYA_SIGNAL_SEQUENCE,
                               PCR_PRODUCT,
                               IDNA,
                               D_CLUSTER,
                               POLYA_SITE,
                               ISRE,
                               MONOCISTRONIC_TRANSCRIPT,
                               DHU_LOOP,
                               THREE_PRIME_D_HEPTAMER,
                               TASIRNA_PRIMARY_TRANSCRIPT,
                               ZINC_FINGER_BINDING_SITE,
                               REPEAT_COMPONENT,
                               LINEAR_DOUBLE_STRANDED_RNA_CHROMOSOME,
                               URIDINE_FIVE_OXYACETIC_ACID,
                               C_TO_G_TRANSVERSION,
                               ASPARTIC_ACID_TRNA_PRIMARY_TRANSCRIPT,
                               CTCF_BINDING_SITE,
                               ANTICODON_LOOP,
                               FOUR_THIOURIDINE,
                               POLYMORPHIC_PSEUDOGENE,
                               R_FIVE_PRIME_LTR_REGION,
                               POLYA_PRIMED_CDNA_CLONE,
                               SPLICE_SITE,
                               DMV2_MOTIF,
                               Z2_REGION,
                               FIVE_PRIME_CODING_EXON_CODING_REGION,
                               FIVE_CARBAMOYLMETHYLURIDINE,
                               FIVE_PRIME_CLIP,
                               CRYPTIC_GENE,
                               SEVEN_METHYLGUANINE,
                               RETROGENE,
                               H3K27_METHYLATION_SITE,
                               NC_PRIMARY_TRANSCRIPT,
                               GROUP_IIB_INTRON,
                               MODIFIED_L_VALINE,
                               H3K79_METHYLATION_SITE,
                               PLASMID_GENE,
                               GENE_SILENCED_BY_RNA_INTERFERENCE,
                               CIS_REGULATORY_FRAMESHIFT_ELEMENT,
                               REGULATORY_REGION,
                               THREE_PRIME_D_RECOMBINATION_SIGNAL_SEQUENCE,
                               MINUS_12_SIGNAL,
                               POSITIVE_SENSE_SSRNA_VIRAL_SEQUENCE,
                               RECOMBINATION_ENHANCER,
                               RNASE_MRP_RNA,
                               GENE_WITH_MRNA_WITH_FRAMESHIFT,
                               ASX_TURN_RIGHT_HANDED_TYPE_TWO,
                               ENHANCERRNA,
                               CHLOROPLAST_CHROMOSOME,
                               FORKHEAD_MOTIF,
                               ATTCTN_SITE,
                               PHAGE_RNA_POLYMERASE_PROMOTER,
                               PSEUDOKNOT,
                               FIVE_METHYLDIHYDROURIDINE,
                               ATTI_SITE,
                               MODIFIED_CYTIDINE,
                               INDEL,
                               D_GENE_SEGMENT,
                               CRYPTIC_SPLICE_SITE,
                               ENGINEERED_REGION,
                               SYMBIOSIS_ISLAND,
                               MIRNA_GENE,
                               ALANYL_TRNA,
                               DYE_TERMINATOR_READ,
                               FIVE_TWO_PRIME_O_DIMETHYLURIDINE,
                               FIVE_METHYLAMINOMETHYL_TWO_THIOURIDINE,
                               TAG,
                               COMPOSITIONALLY_BIASED_REGION_OF_PEPTIDE,
                               RRNA_CLEAVAGE_RNA,
                               POLYPEPTIDE_IRON_ION_CONTACT_SITE,
                               INTERNAL_GUIDE_SEQUENCE,
                               J_C_CLUSTER,
                               ENGINEERED_FOREIGN_REPETITIVE_ELEMENT,
                               ASPARTYL_TRNA,
                               FRAGMENT_ASSEMBLY,
                               CONTIG_COLLECTION,
                               MNP,
                               GROUP_III_INTRON,
                               CLONE_END,
                               RNAPOL_III_PROMOTER_TYPE_3,
                               H3K36_ACETYLATION_SITE,
                               PSEUDOURIDYLATION_GUIDE_SNORNA,
                               PROTEIN_HMM_MATCH,
                               LYSIDINE,
                               CANONICAL_THREE_PRIME_SPLICE_SITE,
                               ATTL_SITE,
                               SIRNA,
                               ALTERNATIVELY_SPLICED_TRANSCRIPT,
                               MODIFIED_L_TRYPTOPHAN,
                               SL9_ACCEPTOR_SITE,
                               REPLICON,
                               MODIFIED_L_TYROSINE,
                               DATABANK_ENTRY,
                               LEFT_HANDED_PEPTIDE_HELIX,
                               CLONE_INSERT_START,
                               MEMBRANE_STRUCTURE,
                               SNORNA,
                               REARRANGEMENT_BREAKPOINT,
                               C_TO_A_TRANSVERSION,
                               H4K91_ACETYLATION_SITE,
                               MODIFIED_L_PROLINE,
                               DECAYED_EXON,
                               EXTENDED_INTRONIC_SPLICE_REGION,
                               V_DJ_C_CLUSTER,
                               GENE_WITH_DICISTRONIC_PRIMARY_TRANSCRIPT,
                               TOPOLOGICALLY_DEFINED_REGION,
                               U5_SNRNA,
                               N2_METHYLGUANOSINE,
                               PROTEIN_MATCH,
                               PLUS_1_TRANSLATIONAL_FRAMESHIFT,
                               PIP_BOX,
                               RRNA_SMALL_SUBUNIT_PRIMARY_TRANSCRIPT,
                               METHYLATED_ADENINE,
                               DUPLICATION,
                               RRNA_5S,
                               H3K4_TRIMETHYLATION,
                               INTERIOR_CODING_EXON,
                               CDS_REGION,
                               RETRON,
                               U12_INTRON,
                               RIGHT_HANDED_PEPTIDE_HELIX,
                               RNAPOL_II_CORE_PROMOTER,
                               SINE_ELEMENT,
                               EXTENDED_CIS_SPLICE_SITE,
                               CRYPTIC_PROPHAGE,
                               METHYLATED_CYTOSINE,
                               STREP_MOTIF,
                               UNTRANSLATED_REGION_POLYCISTRONIC_MRNA,
                               EXTRACHROMOSOMAL_MOBILE_GENETIC_ELEMENT,
                               UTR_REGION,
                               HYDROXYWYBUTOSINE,
                               MODIFIED_GLYCINE,
                               MITOCHONDRIAL_CONTIG,
                               PARALLEL_BETA_STRAND,
                               RIBOSWITCH,
                               NON_PROCESSED_PSEUDOGENE,
                               ENGINEERED_FOREIGN_REGION,
                               DIRECT_TANDEM_DUPLICATION,
                               STOP_CODON_REDEFINED_AS_SELENOCYSTEINE,
                               MICROSATELLITE,
                               TELOMERIC_TRANSCRIPT,
                               T_TO_G_TRANSVERSION,
                               Y_PRIME_ELEMENT,
                               MIRNA_ANTIGUIDE,
                               MINICIRCLE,
                               SL12_ACCEPTOR_SITE,
                               TEMPLATE_REGION,
                               DNA_APTAMER,
                               MITOCHONDRIAL_DNA_READ,
                               POSSIBLE_BASE_CALL_ERROR,
                               PEPTIDE_HELIX,
                               FIVE_METHOXYCARBONYLMETHYLURIDINE,
                               THREE_PRIME_TERMINAL_INVERTED_REPEAT,
                               MODIFIED_L_LYSINE,
                               TCS_ELEMENT,
                               ANTISENSE_LNCRNA,
                               U3_FIVE_PRIME_LTR_REGION,
                               PATERNAL_UNIPARENTAL_DISOMY,
                               PURINE_TRANSITION,
                               TRANSPOSABLE_ELEMENT_GENE,
                               MRNA_REGION,
                               CLEAVED_INITIATOR_METHIONINE,
                               RRNA_16S,
                               EPOXYQUEUOSINE,
                               MIRNA_PRIMARY_TRANSCRIPT,
                               INTEGRON,
                               THREE_METHYLCYTIDINE,
                               EDITING_BLOCK,
                               CHROMOSOME_BAND,
                               CODON,
                               ASX_TURN_LEFT_HANDED_TYPE_ONE,
                               METHYLATION_GUIDE_SNORNA,
                               TRANSCRIPTION_END_SITE,
                               SL1_ACCEPTOR_SITE,
                               U1_SNRNA,
                               CONSTITUTIVE_PROMOTER,
                               RECOMBINATIONALLY_INVERTED_GENE,
                               DCE_SIII,
                               RRNA_25S,
                               MONOCISTRONIC_PRIMARY_TRANSCRIPT,
                               G_QUARTET,
                               LOXP_SITE,
                               FIVE_CARBOXYMETHYLAMINOMETHYLURIDINE,
                               GENOMICALLY_CONTAMINATED_CDNA_CLONE,
                               ISOLEUCINE_TRNA_PRIMARY_TRANSCRIPT,
                               RNAPOL_III_PROMOTER_TYPE_2,
                               MODIFIED_L_ASPARTIC_ACID,
                               PROMOTER_TARGETING_SEQUENCE,
                               TWO_METHYLTHIO_N6_ISOPENTENYLADENOSINE,
                               LINEAR_SINGLE_STRANDED_DNA_CHROMOSOME,
                               Y_REGION,
                               SPLICE_ENHANCER,
                               SMALL_REGULATORY_NCRNA,
                               PIRNA,
                               RDNA_INTERGENIC_SPACER_ELEMENT,
                               DNA_TRANSPOSON,
                               CODING_REGION_OF_EXON,
                               ONE_METHYLINOSINE,
                               DJ_C_CLUSTER,
                               THREE_PRIME_CODING_EXON,
                               PYRIMIDINE_TO_PURINE_TRANSVERSION,
                               TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT,
                               METHYLATION_GUIDE_SNORNA_PRIMARY_TRANSCRIPT,
                               DICISTRONIC_MRNA,
                               REAGENT,
                               PHENYLALANINE_TRNA_PRIMARY_TRANSCRIPT,
                               HISTONE_ACETYLATION_SITE,
                               NON_LTR_RETROTRANSPOSON,
                               ALTERNATE_SEQUENCE_SITE,
                               N2_N2_7_TRIMETHYLGUANOSINE,
                               FIVE_METHYLURIDINE,
                               TSS,
                               FIVE_METHYLAMINOMETHYL_TWO_SELENOURIDINE,
                               ANTISENSE_PRIMARY_TRANSCRIPT,
                               TETRANUCLEOTIDE_REPEAT_MICROSATELLITE_FEATURE,
                               RESTRICTION_ENZYME_FIVE_PRIME_SINGLE_STRAND_OVERHANG,
                               CANONICAL_FIVE_PRIME_SPLICE_SITE,
                               ENHANCER,
                               POLYPEPTIDE_CONSERVED_REGION,
                               TRANSLATION_REGULATORY_REGION,
                               ENGINEERED_FUSION_GENE,
                               GENE_WITH_POLYCISTRONIC_TRANSCRIPT,
                               PHAGE_SEQUENCE,
                               MTE,
                               CPG_ISLAND,
                               POLYPEPTIDE_VARIATION_SITE,
                               KEN_BOX,
                               FIVE_CARBOXYMETHYLAMINOMETHYL_TWO_THIOURIDINE,
                               SL5_ACCEPTOR_SITE,
                               POINT_CENTROMERE,
                               POLYCISTRONIC_PRIMARY_TRANSCRIPT,
                               EXPERIMENTALLY_DEFINED_BINDING_REGION,
                               A_BOX,
                               DIRECT_REPEAT,
                               MATCH_PART,
                               ARCHAEAL_INTRON,
                               R_THREE_PRIME_LTR_REGION,
                               BETA_TURN,
                               DNASEI_HYPERSENSITIVE_SITE,
                               METHIONINE_TRNA_PRIMARY_TRANSCRIPT,
                               PROLINE,
                               STS,
                               INTRAMEMBRANE_POLYPEPTIDE_REGION,
                               YAC_END,
                               V_VDJ_CLUSTER,
                               OLIGO_U_TAIL,
                               JUNCTION,
                               INVERSION_BREAKPOINT,
                               MITOCHONDRIAL_TARGETING_SIGNAL,
                               E_BOX_MOTIF,
                               GENE_SUBARRAY,
                               PROTEIN_PROTEIN_CONTACT,
                               PHEROMONE_RESPONSE_ELEMENT,
                               TRANSLATIONAL_FRAMESHIFT,
                               TERMINAL_INVERTED_REPEAT,
                               A_BOX_TYPE_1,
                               BETA_TURN_RIGHT_HANDED_TYPE_ONE,
                               WYBUTOSINE,
                               RECODED_MRNA,
                               SPLICE_JUNCTION,
                               PROPEPTIDE_CLEAVAGE_SITE,
                               LIPOPROTEIN_SIGNAL_PEPTIDE,
                               VIRTUAL_SEQUENCE,
                               POLYADENYLATED_MRNA,
                               ARIA,
                               AUTOCATALYTICALLY_SPLICED_INTRON,
                               FLEX_ELEMENT,
                               UPSTREAM_AUG_CODON,
                               J_GENE_SEGMENT,
                               AP_1_BINDING_SITE,
                               MORPHOLINO_OLIGO,
                               PROMOTER_TRAP_CONSTRUCT,
                               I_MOTIF,
                               N2_7_DIMETHYLGUANOSINE,
                               RESTRICTION_ENZYME_BINDING_SITE,
                               FIVE_PRIME_TERMINAL_INVERTED_REPEAT,
                               EPITOPE,
                               SRP_RNA,
                               SERYL_TRNA,
                               THREONINE,
                               LIGAND_BINDING_SITE,
                               RECOMBINATION_SIGNAL_SEQUENCE,
                               HISTONE_ACYLATION_REGION,
                               NUCLEOMORPHIC_CHROMOSOME,
                               TRNA_REGION,
                               THREE_PRIME_STEM_LOOP_STRUCTURE,
                               RNA_CHROMOSOME,
                               YAC,
                               HISTONE_METHYLATION_SITE,
                               INOSINE,
                               GNA_OLIGO,
                               RECODED_CODON,
                               ANTI_ARRET,
                               BIOLOGICAL_REGION,
                               PROMOTER_FLANKING_REGION,
                               GROUP_I_INTRON,
                               NON_CANONICAL_START_CODON,
                               PARALOGOUS_REGION,
                               TRANSCRIPTION_PAUSE_SITE,
                               TERMINATOR_OF_TYPE_2_RNAPOL_III_PROMOTER,
                               GENE_WITH_STOP_CODON_READ_THROUGH,
                               INSERTION_SITE,
                               B_BOX,
                               PAIRED_END_FRAGMENT,
                               SERINE_THREONINE_TURN,
                               HISTONE_3_ACETYLATION_SITE,
                               NUCLEOTIDE_BINDING_SITE,
                               N2_7_2PRIRME_O_TRIMETHYLGUANOSINE,
                               GAMMA_TURN_INVERSE,
                               GC_RICH_PROMOTER_REGION,
                               DRE,
                               INTERNAL_SHINE_DALGARNO_SEQUENCE,
                               AFLP_FRAGMENT,
                               SL10_ACCEPTOR_SITE,
                               TWO_THIOURIDINE,
                               PYRROLYSINE_TRNA_PRIMARY_TRANSCRIPT,
                               _5_FORMYLCYTOSINE,
                               INTRONIC_LNCRNA,
                               MRNA_READ,
                               NON_CANONICAL_FIVE_PRIME_SPLICE_SITE,
                               H4K16_ACYLATION_SITE,
                               CODING_EXON,
                               INSERTION_SEQUENCE,
                               PROTEASE_SITE,
                               RRNA_PRIMARY_TRANSCRIPT_REGION,
                               MINI_EXON_DONOR_RNA,
                               RECOMBINATION_FEATURE,
                               SERINE,
                               REGULON,
                               BASE,
                               TR_BOX,
                               BETA_BULGE,
                               VALIDATED_CDNA_CLONE,
                               J_NONAMER,
                               DELETION_BREAKPOINT,
                               RIBOZYME,
                               UTR,
                               GLUTAMINE_TRNA_PRIMARY_TRANSCRIPT,
                               ONE_TWO_PRIME_O_DIMETHYLADENOSINE,
                               ANTISENSE_RNA,
                               N4_N4_2_PRIME_O_TRIMETHYLCYTIDINE,
                               H_ACA_BOX_SNORNA_PRIMARY_TRANSCRIPT,
                               INTERGENIC_REGION,
                               U6ATAC_SNRNA,
                               H3K79_MONOMETHYLATION_SITE,
                               CDS_SUPPORTED_BY_DOMAIN_MATCH_DATA,
                               MIR_ENCODING_TRNA_PRIMARY_TRANSCRIPT,
                               VAULT_RNA,
                               THREE_PRIME_NONCODING_EXON,
                               NESTED_TANDEM_REPEAT,
                               CLONE,
                               R_GNA_OLIGO,
                               METHYLATED_DNA_BASE_FEATURE,
                               ASPARTIC_ACID,
                               NEGATIVELY_AUTOREGULATED_GENE,
                               TWO_METHYLTHIO_N6_CIS_HYDROXYISOPENTENYL_ADENOSINE,
                               C_GENE_SEGMENT,
                               DART_MARKER,
                               MINOR_TSS,
                               RECOMBINATIONALLY_REARRANGED_VERTEBRATE_IMMUNE_SYSTEM_GENE,
                               H2BK20_ACETYLATION_SITE,
                               C_TO_T_TRANSITION_AT_PCPG_SITE,
                               PSEUDOGENE,
                               PATHOGENIC_ISLAND,
                               CHROMOPLAST_CHROMOSOME,
                               EXTERNAL_TRANSCRIBED_SPACER_REGION,
                               RESTRICTION_ENZYME_REGION,
                               FIVE_PRIME_D_RECOMBINATION_SIGNAL_SEQUENCE,
                               STOP_CODON_READ_THROUGH,
                               ADAPTIVE_ISLAND,
                               H2AK5_ACETYLATION_SITE,
                               MIRTRON,
                               H3K36_METHYLATION_SITE,
                               T7_RNA_POLYMERASE_PROMOTER,
                               _2A_SELF_CLEAVING_PEPTIDE_REGION,
                               _2KB_DOWNSTREAM_VARIANT,
                               _2KB_UPSTREAM_VARIANT,
                               _3D_POLYPEPTIDE_STRUCTURE_VARIANT,
                               _3_PRIME_UTR_ELONGATION,
                               _3_PRIME_UTR_EXON_VARIANT,
                               _3_PRIME_UTR_INTRON_VARIANT,
                               _3_PRIME_UTR_TRUNCATION,
                               _3_PRIME_UTR_VARIANT,
                               _500B_DOWNSTREAM_VARIANT,
                               _5KB_DOWNSTREAM_VARIANT,
                               _5KB_UPSTREAM_VARIANT,
                               _5S_SINE_RETROTRANSPOSON,
                               _5_PRIME_UTR_ELONGATION,
                               _5_PRIME_UTR_EXON_VARIANT,
                               _5_PRIME_UTR_INTRON_VARIANT,
                               _5_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT,
                               _5_PRIME_UTR_PREMATURE_START_CODON_LOSS_VARIANT,
                               _5_PRIME_UTR_PREMATURE_START_CODON_VARIANT,
                               _5_PRIME_UTR_TRUNCATION,
                               _5_PRIME_UTR_UORF_FRAMESHIFT_VARIANT,
                               _5_PRIME_UTR_UORF_STOP_CODON_GAIN_VARIANT,
                               _5_PRIME_UTR_UORF_STOP_CODON_LOSS_VARIANT,
                               _5_PRIME_UTR_UORF_STOP_CODON_VARIANT,
                               _5_PRIME_UTR_UORF_VARIANT,
                               _5_PRIME_UTR_VARIANT,
                               _7SL_SINE_RETROTRANSPOSON,
                               AUG_INITIATED_UORF,
                               ALU_DELETION,
                               ALU_INSERTION,
                               BEL_PAO_LTR_RETROTRANSPOSON,
                               C_D_H_ACA_BOX_SCARNA,
                               C_D_H_ACA_BOX_SCARNA_GENE,
                               CAAX_BOX,
                               CACTA_TIR_TRANSPOSON,
                               CARG_BOX,
                               CDS_EXTENSION,
                               CDS_FIVE_PRIME_EXTENSION,
                               CDS_SUPPORTED_BY_PEPTIDE_SPECTRUM_MATCH,
                               CDS_THREE_PRIME_EXTENSION,
                               C_D_BOX_SCARNA,
                               C_D_BOX_SCARNA_GENE,
                               C_D_BOX_SNORNA_ENCODING,
                               C_D_BOX_SNORNA_GENE,
                               COPIA_LTR_RETROTRANSPOSON,
                               CRICK_STRAND,
                               CRYPTON_YR_TRANSPOSON,
                               DIRS_YR_RETROTRANSPOSON,
                               DNA,
                               DNA_LOOP,
                               DNA_LOOP_ANCHOR,
                               DNASEI_HYPERSENSITIVE_SITE_CURRENT_SPELLING,
                               ENDOGENOUS_RETROVIRUS_LTR_RETROTRANSPOSON,
                               FRT_FLANKED,
                               GNA,
                               GYPSY_LTR_RETROTRANSPOSON,
                               H2AZK11_ACETYLATION_SITE,
                               H2AZK13_ACETYLATION_SITE,
                               H2AZK15_ACETYLATION_SITE,
                               H2AZK4_ACETYLATION_SITE,
                               H2AZK7_ACETYLATION_SITE,
                               H2BK5_ACETYLATION_SITE,
                               H3K23_ACETYLATION_SITE,
                               H3K4_ACETYLATION_SITE,
                               H4K16_ACETYLATION_SITE,
                               H4K5_ACETYLATION_SITE,
                               H4K8_ACETYLATION_SITE,
                               HERV_DELETION,
                               HERV_INSERTION,
                               H_ACA_BOX_SCARNA,
                               H_ACA_BOX_SCARNA_GENE,
                               H_ACA_BOX_SNORNA_ENCODING,
                               H_ACA_BOX_SNORNA_GENE,
                               IG_C_GENE,
                               IG_C_PSEUDOGENE,
                               IG_D_GENE,
                               IG_J_GENE,
                               IG_J_PSEUDOGENE,
                               IG_V_GENE,
                               IG_V_PSEUDOGENE,
                               I_LINE_RETROTRANSPOSON,
                               JOCKEY_LINE_RETROTRANSPOSON,
                               L1_LINE_RETROTRANSPOSON,
                               LARD,
                               LINE1_DELETION,
                               LINE1_INSERTION,
                               LNA,
                               LOZ1_RESPONSE_ELEMENT,
                               MNV_ARTIFACT,
                               MAT2P,
                               MAT3M,
                               MERLIN_TIR_TRANSPOSON,
                               MUTATOR_TIR_TRANSPOSON,
                               NMD_ESCAPING_VARIANT,
                               NMD_POLYMORPHIC_PSEUDOGENE_TRANSCRIPT,
                               NMD_TRANSCRIPT,
                               NMD_TRANSCRIPT_VARIANT,
                               NMD_TRIGGERING_VARIANT,
                               NSD_TRANSCRIPT,
                               NGARO_YR_RETROTRANSPOSON,
                               PIF_HARBINGER_TIR_TRANSPOSON,
                               PNA,
                               P_TIR_TRANSPOSON,
                               PENELOPE_RETROTRANSPOSON,
                               PHO7_BINDING_SITE,
                               R2_LINE_RETROTRANSPOSON,
                               RNA,
                               RNA_7SK,
                               RNA_7SK_GENE,
                               RNA_THERMOMETER,
                               RTE_LINE_RETROTRANSPOSON,
                               R_GNA,
                               RETROVIRUS_LTR_RETROTRANSPOSON,
                               ROBERTSONIAN_FUSION,
                               SHP_BOX,
                               SNV_ARTIFACT,
                               SRP_RNA_ENCODING,
                               SUMO_INTERACTION_MOTIF,
                               SVA_DELETION,
                               SVA_INSERTION,
                               S_GNA,
                               TFBS_ABLATION,
                               TFBS_AMPLIFICATION,
                               TFBS_FUSION,
                               TFBS_TRANSLOCATION,
                               TFRS_COLLECTION,
                               TFRS_MODULE,
                               TF_BINDING_SITE_VARIANT,
                               TNA,
                               TRIM,
                               TR_C_GENE,
                               TR_D_GENE,
                               TR_J_GENE,
                               TR_J_PSEUDOGENE,
                               TR_V_GENE,
                               TR_V_PSEUDOGENE,
                               T_CELL_RECEPTOR_GENE,
                               T_CELL_RECEPTOR_PSEUDOGENE,
                               TC1_MARINER_TIR_TRANSPOSON,
                               TRANSIB_TIR_TRANSPOSON,
                               U14_SNORNA_GENE,
                               U3_SNORNA_GENE,
                               U7_SNRNA,
                               UNAAAC_MOTIF,
                               UTR_VARIANT,
                               VIPER_YR_RETROTRANSPOSON,
                               WATSON_STRAND,
                               YR_RETROTRANSPOSON,
                               Y_RNA_GENE,
                               ZAS1_RECOGNITION_MOTIF,
                               ABOLISHED_TRANSCRIPTION_LEVEL_VARIANT,
                               ABSENT_GENE_PRODUCT,
                               ACCESSIBLE_DNA_REGION,
                               ALLELE,
                               ALLELIC_FREQUENCY,
                               ALLELIC_PSEUDOGENE,
                               ALLELIC_PSEUDOGENIC_RRNA,
                               ALLELIC_PSEUDOGENIC_TRNA,
                               ALLELICALLY_EXCLUDED,
                               ALLOPOLYPLOID,
                               ALTERATION_ATTRIBUTE,
                               ALTERED_GENE_PRODUCT_LEVEL,
                               ALTERED_GENE_PRODUCT_SEQUENCE,
                               ALTERNATIVELY_SPLICED,
                               AMBER_STOP_CODON,
                               AMINO_ACID_DELETION,
                               AMINO_ACID_INSERTION,
                               AMINO_ACID_SUBSTITUTION,
                               ANEUPLOID,
                               ANEUPLOID_CHROMOSOME,
                               ANNOTATION_DIRECTED_IMPROVED_DRAFT,
                               ANTISENSE,
                               ANTISENSE_LNCRNA_GENE,
                               APICOPLAST_SEQUENCE,
                               ASSORTMENT_DERIVED_ANEUPLOID,
                               ASSORTMENT_DERIVED_DEFICIENCY,
                               ASSORTMENT_DERIVED_DEFICIENCY_PLUS_DUPLICATION,
                               ASSORTMENT_DERIVED_DUPLICATION,
                               ASSORTMENT_DERIVED_VARIATION,
                               AUTOPOLYPLOID,
                               AUTOREGULATED,
                               AUTOSYNAPTIC_CHROMOSOME,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA54_ELEMENT,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA_70_ELEMENT,
                               BACTERIAL_RNAPOL_PROMOTER_SIGMA_ECF_ELEMENT,
                               BENIGN_VARIANT,
                               BIDIRECTIONAL_GENE_FUSION,
                               BIDIRECTIONAL_PROMOTER_LNCRNA,
                               BIDIRECTIONAL_PROMOTER_LNCRNA_GENE,
                               BIOSYNTHETIC_GENE_CLUSTER,
                               BIPARTITE_DUPLICATION,
                               BIPARTITE_INVERSION,
                               BOUND_BY_FACTOR,
                               BOUND_BY_NUCLEIC_ACID,
                               BOUND_BY_PROTEIN,
                               CDNA,
                               CAPPED,
                               CASSETTE_ARRAY_MEMBER,
                               CELLULAR_MISLOCALIZATION_VARIANT,
                               CHLOROPLAST_DNA,
                               CHLOROPLAST_SEQUENCE,
                               CHROMATIN_REGULATORY_REGION,
                               CHROMOPLAST_SEQUENCE,
                               CHROMOSOMAL_DELETION,
                               CHROMOSOMAL_DUPLICATION,
                               CHROMOSOMAL_INVERSION,
                               CHROMOSOMAL_TRANSLOCATION,
                               CHROMOSOMAL_TRANSPOSITION,
                               CHROMOSOMAL_VARIATION_ATTRIBUTE,
                               CHROMOSOMALLY_ABERRANT_GENOME,
                               CHROMOSOME_FISSION,
                               CHROMOSOME_NUMBER_VARIATION,
                               CHROMOSOME_STRUCTURE_VARIATION,
                               CHROMOSOME_VARIATION,
                               CIRCULAR,
                               CIRCULAR_MRNA,
                               CIRCULAR_NCRNA,
                               CIRCULAR_PLASMID,
                               CIS_REGULATORY_MODULE,
                               CODING_SEQUENCE_VARIANT,
                               CODING_TRANSCRIPT_INTRON_VARIANT,
                               CODING_TRANSCRIPT_VARIANT,
                               CODING_TRANSCRIPT_WITH_RETAINED_INTRON,
                               CODING_VARIANT_QUALITY,
                               CODON_REDEFINED,
                               COMMON_FRAGILE_SITE,
                               COMMON_VARIANT,
                               COMPENSATORY_TRANSCRIPT_SECONDARY_STRUCTURE_VARIANT,
                               COMPLEX_3D_STRUCTURAL_VARIANT,
                               COMPLEX_CHANGE_OF_TRANSLATIONAL_PRODUCT_VARIANT,
                               COMPLEX_CHROMOSOMAL_REARRANGEMENT,
                               COMPLEX_OPERON,
                               COMPLEX_REGULON,
                               COMPLEX_TRANSCRIPT_VARIANT,
                               COMPOUND_CHROMOSOME,
                               COMPOUND_CHROMOSOME_ARM,
                               CONFORMATIONAL_CHANGE_VARIANT,
                               CONSENSUS,
                               CONSERVATIVE_AMINO_ACID_SUBSTITUTION,
                               CONSERVATIVE_INFRAME_DELETION,
                               CONSERVATIVE_INFRAME_INSERTION,
                               CONSERVATIVE_MISSENSE_VARIANT,
                               CONSERVED,
                               CONSERVED_INTERGENIC_VARIANT,
                               CONSERVED_INTRON_VARIANT,
                               COPY_NUMBER_CHANGE,
                               COPY_NUMBER_DECREASE,
                               COPY_NUMBER_INCREASE,
                               CORE_EUKARYOTIC_PROMOTER_ELEMENT,
                               CORE_PROKARYOTIC_PROMOTER_ELEMENT,
                               CORE_VIRAL_PROMOTER_ELEMENT,
                               CRYPTIC,
                               CRYPTIC_PROMOTER,
                               CRYPTIC_SPLICE_ACCEPTOR,
                               CRYPTIC_SPLICE_DONOR,
                               CRYPTIC_SPLICE_SITE_VARIANT,
                               CYANELLE_SEQUENCE,
                               CYCLIC_TRANSLOCATION,
                               CYTOSOLIC_16S_RRNA,
                               CYTOSOLIC_18S_RRNA,
                               CYTOSOLIC_23S_RRNA,
                               CYTOSOLIC_25S_RRNA,
                               CYTOSOLIC_28S_RRNA,
                               CYTOSOLIC_2S_RRNA,
                               CYTOSOLIC_5S_RRNA,
                               CYTOSOLIC_5_8S_RRNA,
                               CYTOSOLIC_LSU_RRNA,
                               CYTOSOLIC_LSU_RRNA_GENE,
                               CYTOSOLIC_SSU_RRNA,
                               CYTOSOLIC_SSU_RRNA_GENE,
                               CYTOSOLIC_RRNA,
                               CYTOSOLIC_RRNA_16S_GENE,
                               CYTOSOLIC_RRNA_18S_GENE,
                               CYTOSOLIC_RRNA_23S_GENE,
                               CYTOSOLIC_RRNA_25S_GENE,
                               CYTOSOLIC_RRNA_28S_GENE,
                               CYTOSOLIC_RRNA_2S_GENE,
                               CYTOSOLIC_RRNA_5S_GENE,
                               CYTOSOLIC_RRNA_5_8S_GENE,
                               CYTOSOLIC_RRNA_GENE,
                               DE_NOVO_VARIANT,
                               DECREASED_ENZYME_ACTIVITY_VARIANT,
                               DECREASED_GENE_PRODUCT_LEVEL,
                               DECREASED_POLYADENYLATION_VARIANT,
                               DECREASED_POLYPEPTIDE_LEVEL,
                               DECREASED_TRANSCRIPT_LEVEL_VARIANT,
                               DECREASED_TRANSCRIPT_STABILITY_VARIANT,
                               DECREASED_TRANSCRIPTION_RATE_VARIANT,
                               DECREASED_TRANSLATIONAL_PRODUCT_LEVEL,
                               DEFICIENT_INTERCHROMOSOMAL_TRANSPOSITION,
                               DEFICIENT_INTRACHROMOSOMAL_TRANSPOSITION,
                               DEFICIENT_INVERSION,
                               DEFICIENT_TRANSLOCATION,
                               DELETION_ARTIFACT,
                               DELINS,
                               DESIGNED_SEQUENCE,
                               DEXTROSYNAPTIC_CHROMOSOME,
                               DICISTRONIC,
                               DIPLOTYPE,
                               DIRECT,
                               DIRECTION_ATTRIBUTE,
                               DISEASE_ASSOCIATED_VARIANT,
                               DISEASE_CAUSING_VARIANT,
                               DISRUPTIVE_INFRAME_DELETION,
                               DISRUPTIVE_INFRAME_INSERTION,
                               DOMINANT_NEGATIVE_VARIANT,
                               DOUBLE,
                               DOUBLE_STRANDED_CDNA,
                               DOWNSTREAM_GENE_VARIANT,
                               DOWNSTREAM_INTERGENIC_FUSION,
                               DOWNSTREAM_TRANSCRIPT_VARIANT,
                               DUPLICATION_ARTIFACT,
                               DUPLICATION_ATTRIBUTE,
                               EARLY_ORIGIN_OF_REPLICATION,
                               EDITED,
                               EDITING_VARIANT,
                               ELONGATED_IN_FRAME_POLYPEPTIDE_C_TERMINAL,
                               ELONGATED_IN_FRAME_POLYPEPTIDE_N_TERMINAL_ELONGATION,
                               ELONGATED_OUT_OF_FRAME_POLYPEPTIDE_C_TERMINAL,
                               ELONGATED_OUT_OF_FRAME_POLYPEPTIDE_N_TERMINAL,
                               ELONGATED_POLYPEPTIDE,
                               ELONGATED_POLYPEPTIDE_C_TERMINAL,
                               ELONGATED_POLYPEPTIDE_N_TERMINAL,
                               ENCODES_1_POLYPEPTIDE,
                               ENCODES_ALTERNATE_TRANSCRIPTION_START_SITES,
                               ENCODES_ALTERNATELY_SPLICED_TRANSCRIPTS,
                               ENCODES_DIFFERENT_POLYPEPTIDES_DIFFERENT_STOP,
                               ENCODES_DISJOINT_POLYPEPTIDES,
                               ENCODES_GREATER_THAN_1_POLYPEPTIDE,
                               ENCODES_OVERLAPPING_PEPTIDES,
                               ENCODES_OVERLAPPING_PEPTIDES_DIFFERENT_START,
                               ENCODES_OVERLAPPING_POLYPEPTIDES_DIFFERENT_START_AND_STOP,
                               ENDOGENOUS_RETROVIRAL_SEQUENCE,
                               ENGINEERED,
                               ENHANCER_BLOCKING_ELEMENT,
                               ENZYMATIC,
                               ENZYMATIC_RNA_GENE,
                               EPIGENETICALLY_MODIFIED,
                               EPIGENOMICALLY_MODIFIED_REGION,
                               EUKARYOTIC_PROMOTER,
                               EXEMPLAR,
                               EXON_LOSS_VARIANT,
                               EXON_VARIANT,
                               EXONIC_SPLICE_REGION_VARIANT,
                               EXPERIMENTAL_FEATURE_ATTRIBUTE,
                               EXPERIMENTALLY_DETERMINED,
                               EXTENDED_INTRONIC_SPLICE_REGION_VARIANT,
                               FEATURE_ABLATION,
                               FEATURE_AMPLIFICATION,
                               FEATURE_ATTRIBUTE,
                               FEATURE_ELONGATION,
                               FEATURE_FUSION,
                               FEATURE_TRANSLOCATION,
                               FEATURE_TRUNCATION,
                               FEATURE_VARIANT,
                               FINISHED_GENOME,
                               FIVE_PRIME_UTR_PREMATURE_START_CODON_LOCATION_VARIANT,
                               FIVE_PRIME_DUPLICATED_TRANSCRIPT,
                               FIVE_PRIME_FIVE_PRIME_OVERLAP,
                               FIVE_PRIME_GENE_FUSION,
                               FIVE_PRIME_THREE_PRIME_OVERLAP,
                               FIXED_VARIANT,
                               FLANKED,
                               FLANKING_REPEAT,
                               FLOXED,
                               FOREIGN,
                               FORWARD,
                               FRAGILE_SITE,
                               FRAGMENTARY,
                               FRAME_RESTORING_VARIANT,
                               FRAMESHIFT,
                               FRAMESHIFT_ELONGATION,
                               FRAMESHIFT_TRUNCATION,
                               FRAMESHIFT_UNIDIRECTIONAL_GENE_FUSION,
                               FRAMESHIFT_VARIANT,
                               FRAMESHIFT_VARIANT_NMD_ESCAPING,
                               FRAMESHIFT_VARIANT_NMD_TRIGGERING,
                               FREE,
                               FREE_CHROMOSOME_ARM,
                               FREE_DUPLICATION,
                               FREE_RING_DUPLICATION,
                               FUNCTION_UNCERTAIN_VARIANT,
                               FUNCTIONAL_EFFECT_VARIANT,
                               FUNCTIONAL_GENE_REGION,
                               FUNCTIONALLY_ABNORMAL,
                               FUNCTIONALLY_NORMAL,
                               FUSION,
                               GRNA_ENCODING,
                               GAIN_OF_FUNCTION_VARIANT,
                               GENE_ARRAY_MEMBER,
                               GENE_ATTRIBUTE,
                               GENE_CASSETTE_MEMBER,
                               GENE_FUSION,
                               GENE_SUBARRAY_MEMBER,
                               GENE_TO_GENE_FEATURE,
                               GENE_VARIANT,
                               GENIC_DOWNSTREAM_TRANSCRIPT_VARIANT,
                               GENIC_UPSTREAM_TRANSCRIPT_VARIANT,
                               GENOME,
                               GENOMIC_DNA,
                               GENOMIC_DNA_CONTIG,
                               GENOMICALLY_IMPRINTED,
                               GENOTYPE,
                               GERMLINE_VARIANT,
                               GROUP_IIC_INTRON,
                               HAT_TIR_TRANSPOSON,
                               HAPLOTYPE,
                               HETERO_COMPOUND_CHROMOSOME,
                               HIGH_QUALITY_DRAFT,
                               HISTONE_2AZ_ACETYLATION_SITE,
                               HISTONE_2A_ACETYLATION_SITE,
                               HISTONE_2B_ACETYLATION_SITE,
                               HISTONE_4_ACETYLATION_SITE,
                               HOMO_COMPOUND_CHROMOSOME,
                               HOMOLOGOUS,
                               HOMOLOGOUS_CHROMOSOME_RECOGNITION_AND_PAIRING_LOCUS,
                               HPRNA,
                               HPRNA_GENE,
                               HYPERPLOID,
                               HYPOPLOID,
                               IMMUNOGLOBULIN_GENE,
                               IMMUNOGLOBULIN_PSEUDOGENE,
                               IMPRINTING_CONTROL_REGION,
                               IMPROVED_HIGH_QUALITY_DRAFT,
                               INACTIVE_CATALYTIC_SITE,
                               INACTIVE_LIGAND_BINDING_SITE,
                               INCOMPLETE_TERMINAL_CODON_VARIANT,
                               INCOMPLETE_TRANSCRIPT_3UTR_VARIANT,
                               INCOMPLETE_TRANSCRIPT_5UTR_VARIANT,
                               INCOMPLETE_TRANSCRIPT_CDS,
                               INCOMPLETE_TRANSCRIPT_CODING_SPLICE_VARIANT,
                               INCOMPLETE_TRANSCRIPT_EXONIC_VARIANT,
                               INCOMPLETE_TRANSCRIPT_INTRONIC_VARIANT,
                               INCOMPLETE_TRANSCRIPT_SPLICE_REGION_VARIANT,
                               INCOMPLETE_TRANSCRIPT_VARIANT,
                               INCREASED_ENZYME_ACTIVITY_VARIANT,
                               INCREASED_GENE_PRODUCT_LEVEL,
                               INCREASED_POLYADENYLATION_VARIANT,
                               INCREASED_TRANSCRIPT_LEVEL_VARIANT,
                               INCREASED_TRANSCRIPT_STABILITY_VARIANT,
                               INCREASED_TRANSCRIPTION_RATE_VARIANT,
                               INCREASED_TRANSLATIONAL_PRODUCT_LEVEL,
                               INDEL_ARTIFACT,
                               INDEPENDENTLY_KNOWN,
                               INERT_DNA_SPACER,
                               INFRAME,
                               INFRAME_DELETION,
                               INFRAME_INDEL,
                               INFRAME_INSERTION,
                               INFRAME_UNIDIRECTIONAL_GENE_FUSION,
                               INFRAME_VARIANT,
                               INITIATOR_CODON_VARIANT,
                               INSERTION_ARTIFACT,
                               INSERTION_ATTRIBUTE,
                               INSERTIONAL,
                               INSERTIONAL_DUPLICATION,
                               INSIDE_INTRON,
                               INSIDE_INTRON_ANTIPARALLEL,
                               INSIDE_INTRON_PARALLEL,
                               INTEIN_CONTAINING,
                               INTERCHROMOSOMAL,
                               INTERCHROMOSOMAL_DUPLICATION,
                               INTERCHROMOSOMAL_MUTATION,
                               INTERCHROMOSOMAL_TRANSLOCATION,
                               INTERCHROMOSOMAL_TRANSPOSITION,
                               INTERGENIC_1KB_VARIANT,
                               INTERGENIC_VARIANT,
                               INTERMEDIATE,
                               INTERNAL_FEATURE_ELONGATION,
                               INTRACHROMOSOMAL,
                               INTRACHROMOSOMAL_DUPLICATION,
                               INTRACHROMOSOMAL_MUTATION,
                               INTRACHROMOSOMAL_TRANSLOCATION,
                               INTRACHROMOSOMAL_TRANSPOSITION,
                               INTRAGENIC_VARIANT,
                               INTRON_GAIN_VARIANT,
                               INTRON_VARIANT,
                               INVALIDATED,
                               INVALIDATED_BY_CHIMERIC_CDNA,
                               INVALIDATED_BY_GENOMIC_CONTAMINATION,
                               INVALIDATED_BY_GENOMIC_POLYA_PRIMED_CDNA,
                               INVALIDATED_BY_PARTIAL_PROCESSING,
                               INVERSION_ATTRIBUTE,
                               INVERSION_CUM_TRANSLOCATION,
                               INVERSION_DERIVED_ANEUPLOID_CHROMOSOME,
                               INVERSION_DERIVED_BIPARTITE_DEFICIENCY,
                               INVERSION_DERIVED_BIPARTITE_DUPLICATION,
                               INVERSION_DERIVED_DEFICIENCY_PLUS_ANEUPLOID,
                               INVERSION_DERIVED_DEFICIENCY_PLUS_DUPLICATION,
                               INVERSION_DERIVED_DUPLICATION_PLUS_ANEUPLOID,
                               INVERTED,
                               INVERTED_INSERTIONAL_DUPLICATION,
                               INVERTED_INTERCHROMOSOMAL_TRANSPOSITION,
                               INVERTED_INTRACHROMOSOMAL_TRANSPOSITION,
                               INVERTED_RING_CHROMOSOME,
                               ISOMIR,
                               KINETOPLAST,
                               KNOB,
                               LAEVOSYNAPTIC_CHROMOSOME,
                               LATE_ORIGIN_OF_REPLICATION,
                               LETHAL_VARIANT,
                               LEUCOPLAST_SEQUENCE,
                               LEVEL_OF_TRANSCRIPT_VARIANT,
                               LINEAR,
                               LINEAR_PLASMID,
                               LNCRNA_GENE,
                               LNCRNA_WITH_RETAINED_INTRON,
                               LONG_TERMINAL_REPEAT_TRANSCRIPT,
                               LOOP,
                               LOSS_OF_FUNCTION_VARIANT,
                               LOSS_OF_HETEROZYGOSITY,
                               LOST_POLYPEPTIDE,
                               LOW_COMPLEXITY,
                               MRNA_ATTRIBUTE,
                               MACRONUCLEAR_SEQUENCE,
                               MATERNAL_VARIANT,
                               MATERNALLY_IMPRINTED,
                               MATURE_MIRNA_VARIANT,
                               MATURE_PROTEIN_REGION_OF_CDS,
                               MEIOTIC_RECOMBINATION_REGION,
                               MEMBER_OF_REGULON,
                               METHYLATION_GUIDE_SNORNA_GENE,
                               MIRNA_ENCODING,
                               MICRONUCLEAR_SEQUENCE,
                               MINUS_1_FRAMESHIFT,
                               MINUS_1_FRAMESHIFT_VARIANT,
                               MINUS_1_TRANSLATIONAL_FRAMESHIFT,
                               MINUS_1_TRANSLATIONALLY_FRAMESHIFTED,
                               MINUS_2_FRAMESHIFT,
                               MINUS_2_FRAMESHIFT_VARIANT,
                               MINUS_2_TRANSLATIONAL_FRAMESHIFT,
                               MISSENSE_VARIANT,
                               MITOCHONDRIAL_DNA,
                               MITOCHONDRIAL_D_LOOP,
                               MITOCHONDRIAL_CONTROL_REGION,
                               MITOCHONDRIAL_SEQUENCE,
                               MITOTIC_RECOMBINATION_REGION,
                               MOBILE,
                               MOBILE_ELEMENT_DELETION,
                               MONOCISTRONIC,
                               MORPHOLINO_BACKBONE,
                               MT_LSU_RRNA,
                               MT_LSU_RRNA_GENE,
                               MT_SSU_RRNA,
                               MT_SSU_RRNA_GENE,
                               MT_RRNA,
                               MT_RRNA_GENE,
                               MT_TRNA,
                               MUTATIONAL_HOTSPOT,
                               NATURAL,
                               NEGATIVELY_AUTOREGULATED,
                               NO_SEQUENCE_ALTERATION,
                               NON_AUG_INITIATED_UORF,
                               NON_ALLELIC_HOMOLOGOUS_RECOMBINATION_REGION,
                               NON_CODING_TRANSCRIPT_EXON_VARIANT,
                               NON_CODING_TRANSCRIPT_INTRON_VARIANT,
                               NON_CODING_TRANSCRIPT_SPLICE_REGION_VARIANT,
                               NON_CODING_TRANSCRIPT_VARIANT,
                               NON_COMPLIMENTARY_STEM,
                               NON_CONSERVATIVE_AMINO_ACID_SUBSTITUTION,
                               NON_CONSERVATIVE_MISSENSE_VARIANT,
                               NON_PROTEIN_CODING,
                               NON_SYNONYMOUS,
                               NONCONTIGUOUS_FINISHED,
                               NONSYNONYMOUS_VARIANT,
                               NUCLEAR_SEQUENCE,
                               NUCLEIC_ACID,
                               NUCLEOMORPHIC_SEQUENCE,
                               NUCLEOTIDE_CLEAVAGE_SITE,
                               NULL_MUTATION,
                               OCHRE_STOP_CODON,
                               ONCOGENIC_VARIANT,
                               OPAL_STOP_CODON,
                               OPERON_MEMBER,
                               ORGANELLE_SEQUENCE,
                               ORPHAN,
                               ORTHOLOGOUS,
                               OVERLAPPING,
                               PARACENTRIC,
                               PARACENTRIC_INVERSION,
                               PARALOGOUS,
                               PARTIALLY_CHARACTERIZED_CHROMOSOMAL_MUTATION,
                               PARTIALLY_DUPLICATED_TRANSCRIPT,
                               PATERNAL_VARIANT,
                               PATERNALLY_IMPRINTED,
                               PEDIGREE_SPECIFIC_VARIANT,
                               PEPTIDE_COLLECTION,
                               PEPTIDYL,
                               PERICENTRIC,
                               PERICENTRIC_INVERSION,
                               PIGGYBAC_TIR_TRANSPOSON,
                               PLASMID_LOCATION,
                               PLASTID_LSU_RRNA,
                               PLASTID_LSU_RRNA_GENE,
                               PLASTID_SSU_RRNA,
                               PLASTID_SSU_RRNA_GENE,
                               PLASTID_RRNA,
                               PLASTID_RRNA_GENE,
                               PLASTID_SEQUENCE,
                               PLUS_1_FRAMESHIFT,
                               PLUS_1_FRAMESHIFT_VARIANT,
                               PLUS_1_TRANSLATIONALLY_FRAMESHIFTED,
                               PLUS_2_FRAMESHIFT_VARIANT,
                               PLUS_2_FRAMSHIFT,
                               POLYA_SITE_CLUSTER,
                               POLYADENYLATED,
                               POLYADENYLATION_VARIANT,
                               POLYCISTRONIC,
                               POLYMER_ATTRIBUTE,
                               POLYMORPHIC_PSEUDOGENE_PROCESSED_TRANSCRIPT,
                               POLYMORPHIC_PSEUDOGENE_WITH_RETAINED_INTRON,
                               POLYMORPHIC_SEQUENCE_VARIANT,
                               POLYMORPHIC_VARIANT,
                               POLYPEPTIDE_FUNCTION_VARIANT,
                               POLYPEPTIDE_FUSION,
                               POLYPEPTIDE_GAIN_OF_FUNCTION_VARIANT,
                               POLYPEPTIDE_LOCALIZATION_VARIANT,
                               POLYPEPTIDE_LOSS_OF_FUNCTION_VARIANT,
                               POLYPEPTIDE_PARTIAL_LOSS_OF_FUNCTION,
                               POLYPEPTIDE_POST_TRANSLATIONAL_PROCESSING_VARIANT,
                               POLYPEPTIDE_SEQUENCE_VARIANT,
                               POLYPEPTIDE_TRUNCATION,
                               POLYPLOID,
                               POPULATION_SPECIFIC_VARIANT,
                               POSITIVELY_AUTOREGULATED,
                               POST_TRANSLATIONALLY_REGULATED,
                               POST_TRANSLATIONALLY_REGULATED_BY_PROTEIN_MODIFICATION,
                               POST_TRANSLATIONALLY_REGULATED_BY_PROTEIN_STABILITY,
                               PRECURSOR_RNA,
                               PREDICTED,
                               PREDICTED_BY_AB_INITIO_COMPUTATION,
                               PREDICTED_TRANSCRIPT,
                               PRESENCE_ABSENCE_VARIATION,
                               PROCESSED_PSEUDOGENIC_RRNA,
                               PROCESSED_PSEUDOGENIC_TRNA,
                               PROKARYOTIC_PROMOTER,
                               PROPEPTIDE_REGION_OF_CDS,
                               PROPLASTID_SEQUENCE,
                               PROTEIN_ALTERING_VARIANT,
                               PROTEIN_CODING,
                               PROVIRAL_LOCATION,
                               PSEUDOEXON_VARIANT,
                               PSEUDOGENE_PROCESSED_TRANSCRIPT,
                               PSEUDOGENIC_CDS,
                               PSEUDOGENIC_TRANSCRIPT_WITH_RETAINED_INTRON,
                               PSEUDOURIDYLATION_GUIDE_SNORNA_GENE,
                               PUMILIO_RESPONSE_ELEMENT,
                               PYRROLYSINE_LOSS,
                               QUALITY_VALUE,
                               QUANTITATIVE_VARIANT,
                               RRNA_ENCODING,
                               RANDOM_SEQUENCE,
                               RANGE_EXTENDER_ELEMENT,
                               RARE_AMINO_ACID_VARIANT,
                               RARE_FRAGILE_SITE,
                               RARE_VARIANT,
                               RATE_OF_TRANSCRIPTION_VARIANT,
                               REARRANGED_AT_DNA_LEVEL,
                               REARRANGEMENT_REGION,
                               RECIPROCAL,
                               RECIPROCAL_CHROMOSOMAL_TRANSLOCATION,
                               RECODED,
                               RECODED_BY_TRANSLATIONAL_BYPASS,
                               RECOMBINATIONALLY_REARRANGED,
                               REDUNDANT_INSERTED_STOP_GAINED,
                               REF_MIRNA,
                               REFERENCE_GENOME,
                               REGULATED,
                               REGULATORY_REGION_ABLATION,
                               REGULATORY_REGION_AMPLIFICATION,
                               REGULATORY_REGION_FUSION,
                               REGULATORY_REGION_TRANSLOCATION,
                               REGULATORY_REGION_VARIANT,
                               REPEAT_INSTABILITY_REGION,
                               REPLICATION_START_SITE,
                               RESCUE,
                               RESPONSE_ELEMENT,
                               RETROTRANSPOSED,
                               REVERSE,
                               RIBOZYME_GENE,
                               RIBOZYMIC,
                               RING_CHROMOSOME,
                               SBRNA,
                               SBRNA_GENE,
                               SCRNA_ENCODING,
                               SCARNA,
                               SCARNA_GENE,
                               SCORE,
                               SELECTION_MARKER,
                               SELENOCYSTEINE_GAIN,
                               SELENOCYSTEINE_LOSS,
                               SELF_CLEAVING_RIBOZYME,
                               SENSE_INTRONIC_LNCRNA,
                               SENSE_INTRONIC_LNCRNA_GENE,
                               SENSE_OVERLAP_LNCRNA,
                               SENSE_OVERLAP_LNCRNA_GENE,
                               SEQUENCE_ALTERATION_ARTIFACT,
                               SEQUENCE_ATTRIBUTE,
                               SEQUENCE_COLLECTION,
                               SEQUENCE_COMPARISON,
                               SEQUENCE_LENGTH_ALTERATION,
                               SEQUENCE_LENGTH_VARIANT,
                               SEQUENCE_LOCATION,
                               SEQUENCE_SOURCE,
                               SEQUENCE_VARIANT,
                               SHORT_TANDEM_REPEAT_CHANGE,
                               SHORT_TANDEM_REPEAT_CONTRACTION,
                               SHORT_TANDEM_REPEAT_EXPANSION,
                               SHORT_TANDEM_REPEAT_VARIATION,
                               SIGNAL_PEPTIDE_REGION_OF_CDS,
                               SILENCED,
                               SILENCED_BY_DNA_METHYLATION,
                               SILENCED_BY_DNA_MODIFICATION,
                               SILENCED_BY_RNA_INTERFERENCE,
                               SILENCED_BY_HISTONE_DEACETYLATION,
                               SILENCED_BY_HISTONE_METHYLATION,
                               SILENCED_BY_HISTONE_MODIFICATION,
                               SILENT_MUTATION,
                               SIMPLE_OPERON,
                               SIMPLE_REGULON,
                               SINGLE,
                               SINGLE_STRANDED_CDNA,
                               SISRNA,
                               SNRNA_ENCODING,
                               SNCRNA,
                               SNCRNA_GENE,
                               SNORNA_ENCODING,
                               SOMATIC_VARIANT,
                               SPLICE_ACCEPTOR_VARIANT,
                               SPLICE_ACCEPTOR_VARIANT_NMD_ESCAPING,
                               SPLICE_ACCEPTOR_VARIANT_NMD_TRIGGERING,
                               SPLICE_DONOR_5TH_BASE_VARIANT,
                               SPLICE_DONOR_REGION_VARIANT,
                               SPLICE_DONOR_VARIANT,
                               SPLICE_DONOR_VARIANT_NMD_ESCAPING,
                               SPLICE_DONOR_VARIANT_NMD_TRIGGERING,
                               SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
                               SPLICE_REGION_VARIANT,
                               SPLICE_SITE_VARIANT,
                               SPLICING_VARIANT,
                               SPURIOUS_PROTEIN,
                               STRNA_ENCODING,
                               STANDARD_DRAFT,
                               START_LOST,
                               START_RETAINED_VARIANT,
                               STATUS,
                               STEM,
                               STEM_LOOP_REGION,
                               STOP_GAINED,
                               STOP_GAINED_NMD_ESCAPING,
                               STOP_GAINED_NMD_TRIGGERING,
                               STOP_LOST,
                               STOP_RETAINED_VARIANT,
                               STRAND_ATTRIBUTE,
                               STRUCTURAL_INTERACTION_VARIANT,
                               STRUCTURAL_VARIANT,
                               SUBSTITUTION_ARTIFACT,
                               SUPPORTED_BY_EST_OR_CDNA,
                               SUPPORTED_BY_DOMAIN_MATCH,
                               SUPPORTED_BY_SEQUENCE_SIMILARITY,
                               SYNONYMOUS,
                               SYNONYMOUS_VARIANT,
                               SYNTENIC,
                               SYNTHETIC_SEQUENCE,
                               TRNA_SINE_RETROTRANSPOSON,
                               TRNA_ENCODING,
                               TANDEM,
                               TEB1_RECOGNITION_MOTIF,
                               TELOMERIC_D_LOOP,
                               TERMINAL_REPEAT,
                               TERMINATOR_CODON_VARIANT,
                               THREE_PRIME_DUPLICATED_TRANSCRIPT,
                               THREE_PRIME_FIVE_PRIME_OVERLAP,
                               THREE_PRIME_GENE_FUSION,
                               THREE_PRIME_OVERLAPPING_NCRNA,
                               THREE_PRIME_THREE_PRIME_OVERLAP,
                               TMRNA_ENCODING,
                               TOPOLOGICALLY_ASSOCIATED_DOMAIN,
                               TOPOLOGICALLY_ASSOCIATED_DOMAIN_BOUNDARY,
                               TOPOLOGY_ATTRIBUTE,
                               TRANS_SPLICED,
                               TRANSCRIBED_PROCESSED_PSEUDOGENE,
                               TRANSCRIBED_UNITARY_PSEUDOGENE,
                               TRANSCRIBED_UNPROCESSED_PSEUDOGENE,
                               TRANSCRIPT_ABLATION,
                               TRANSCRIPT_AMPLIFICATION,
                               TRANSCRIPT_ATTRIBUTE,
                               TRANSCRIPT_FUNCTION_VARIANT,
                               TRANSCRIPT_FUSION,
                               TRANSCRIPT_PROCESSING_VARIANT,
                               TRANSCRIPT_REGULATORY_REGION_FUSION,
                               TRANSCRIPT_SECONDARY_STRUCTURE_VARIANT,
                               TRANSCRIPT_STABILITY_VARIANT,
                               TRANSCRIPT_TRANSLOCATION,
                               TRANSCRIPT_VARIANT,
                               TRANSCRIPTION_FACTOR_REGULATORY_SITE,
                               TRANSCRIPTION_TERMINATION_SIGNAL,
                               TRANSCRIPTION_UNIT,
                               TRANSCRIPTION_VARIANT,
                               TRANSCRIPTIONALLY_CONSTITUTIVE,
                               TRANSCRIPTIONALLY_INDUCED,
                               TRANSCRIPTIONALLY_REGULATED,
                               TRANSCRIPTIONALLY_REPRESSED,
                               TRANSGENIC,
                               TRANSIT_PEPTIDE_REGION_OF_CDS,
                               TRANSLATED_PROCESSED_PSEUDOGENE,
                               TRANSLATED_UNPROCESSED_PSEUDOGENE,
                               TRANSLATIONAL_PRODUCT_FUNCTION_VARIANT,
                               TRANSLATIONAL_PRODUCT_LEVEL_VARIANT,
                               TRANSLATIONAL_PRODUCT_STRUCTURE_VARIANT,
                               TRANSLATIONALLY_FRAMESHIFTED,
                               TRANSLATIONALLY_REGULATED,
                               TRANSLOCATION_ELEMENT,
                               TRANSLOCATON_ATTRIBUTE,
                               TRINUCLEOTIDE_REPEAT_EXPANSION,
                               UNCHARACTERIZED_CHROMOSOMAL_MUTATION,
                               UNCONFIRMED_TRANSCRIPT,
                               UNIDIRECTIONAL_GENE_FUSION,
                               UNIDIRECTIONAL_INTERGENIC_FUSION,
                               UNINVERTED_INSERTIONAL_DUPLICATION,
                               UNINVERTED_INTERCHROMOSOMAL_TRANSPOSITION,
                               UNINVERTED_INTRACHROMOSOMAL_TRANSPOSITION,
                               UNIQUE_VARIANT,
                               UNIT_OF_GENE_EXPRESSION,
                               UNITARY_PSEUDOGENIC_RRNA,
                               UNITARY_PSEUDOGENIC_TRNA,
                               UNORIENTED_INSERTIONAL_DUPLICATION,
                               UNORIENTED_INTERCHROMOSOMAL_TRANSPOSITION,
                               UNORIENTED_INTRACHROMOSOMAL_TRANSPOSITION,
                               UNPROCESSED_PSEUDOGENIC_RRNA,
                               UNPROCESSED_PSEUDOGENIC_TRNA,
                               UNSPECIFIED_INDEL,
                               UPSTREAM_GENE_VARIANT,
                               UPSTREAM_INTERGENIC_FUSION,
                               UPSTREAM_TRANSCRIPT_VARIANT,
                               VALIDATED,
                               VARIANT_COLLECTION,
                               VARIANT_FREQUENCY,
                               VARIANT_GENOME,
                               VARIANT_ORIGIN,
                               VARIANT_PHENOTYPE,
                               VARIANT_QUALITY,
                               VAULT_RNA_GENE,
                               VERTEBRATE_IMMUNE_SYSTEM_GENE,
                               VERTEBRATE_IMMUNE_SYSTEM_PSEUDOGENE,
                               VIRAL_PROMOTER,
                               WHOLE_GENOME_SEQUENCE_STATUS,
                               WILD_TYPE,
                               })

GENE = "gene"
EPIGENETICALLY_MODIFIED_GENE = "epigenetically_modified_gene"
GENE_CASSETTE = "gene_cassette"
GENE_WITH_START_CODON_CUG = "gene_with_start_codon_CUG"
GENE_WITH_TRANS_SPLICED_TRANSCRIPT = "gene_with_trans_spliced_transcript"
TRANSLATIONALLY_REGULATED_GENE = "translationally_regulated_gene"
CT_GENE = "ct_gene"
GENE_WITH_STOP_CODON_REDEFINED_AS_SELENOCYSTEINE = "gene_with_stop_codon_redefined_as_selenocysteine"
RECOMBINATIONALLY_INVERTED_GENE = "recombinationally_inverted_gene"
PLASTID_GENE = "plastid_gene"
GENE_WITH_DICISTRONIC_MRNA = "gene_with_dicistronic_mRNA"
FOREIGN_GENE = "foreign_gene"
GENE_WITH_RECODED_MRNA = "gene_with_recoded_mRNA"
CRYPTOGENE = "cryptogene"
KINETOPLAST_GENE = "kinetoplast_gene"
APICOPLAST_GENE = "apicoplast_gene"
CHROMOPLAST_GENE = "chromoplast_gene"
GENE_WITH_STOP_CODON_READ_THROUGH = "gene_with_stop_codon_read_through"
GENE_WITH_TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT = "gene_with_transcript_with_translational_frameshift"
GENE_REARRANGED_AT_DNA_LEVEL = "gene_rearranged_at_DNA_level"
NUCLEOMORPH_GENE = "nucleomorph_gene"
TRANSGENE = "transgene"
FUNCTIONAL_CANDIDATE_GENE = "functional_candidate_gene"
GENE_WITH_DICISTRONIC_PRIMARY_TRANSCRIPT = "gene_with_dicistronic_primary_transcript"
PIRNA_GENE = "piRNA_gene"
ENGINEERED_GENE = "engineered_gene"
POSITIONAL_CANDIDATE_GENE = "positional_candidate_gene"
CANDIDATE_GENE = "candidate_gene"
RNASE_MRP_RNA_GENE = "RNase_MRP_RNA_gene"
GENE_WITH_STOP_CODON_REDEFINED_AS_PYRROLYSINE = "gene_with_stop_codon_redefined_as_pyrrolysine"
PROVIRAL_GENE = "proviral_gene"
SRP_RNA_GENE = "SRP_RNA_gene"
TRNA_GENE = "tRNA_gene"
CRYPTIC_GENE = "cryptic_gene"
FUSION_GENE = "fusion_gene"
ALLELICALLY_EXCLUDED_GENE = "allelically_excluded_gene"
RETROGENE = "retrogene"
GRNA_GENE = "gRNA_gene"
ENGINEERED_FUSION_GENE = "engineered_fusion_gene"
GENE_WITH_POLYCISTRONIC_TRANSCRIPT = "gene_with_polycistronic_transcript"
PLASMID_GENE = "plasmid_gene"
MATERNALLY_IMPRINTED_GENE = "maternally_imprinted_gene"
GENE_SILENCED_BY_RNA_INTERFERENCE = "gene_silenced_by_RNA_interference"
WILD_TYPE_RESCUE_GENE = "wild_type_rescue_gene"
RNASE_P_RNA_GENE = "RNase_P_RNA_gene"
RRNA_GENE = "rRNA_gene"
GENE_SILENCED_BY_DNA_MODIFICATION = "gene_silenced_by_DNA_modification"
ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT_GENE = "engineered_foreign_transposable_element_gene"
SILENCED_GENE = "silenced_gene"
GENE_WITH_POLYADENYLATED_MRNA = "gene_with_polyadenylated_mRNA"
MINICIRCLE_GENE = "minicircle_gene"
PREDICTED_GENE = "predicted_gene"
TELOMERASE_RNA_GENE = "telomerase_RNA_gene"
RESCUE_GENE = "rescue_gene"
PATERNALLY_IMPRINTED_GENE = "paternally_imprinted_gene"
NUCLEAR_GENE = "nuclear_gene"
GENE_WITH_MRNA_WITH_FRAMESHIFT = "gene_with_mRNA_with_frameshift"
SCRNA_GENE = "scRNA_gene"
NEGATIVELY_AUTOREGULATED_GENE = "negatively_autoregulated_gene"
GENE_SILENCED_BY_HISTONE_METHYLATION = "gene_silenced_by_histone_methylation"
GENE_SILENCED_BY_HISTONE_MODIFICATION = "gene_silenced_by_histone_modification"
SNRNA_GENE = "snRNA_gene"
POST_TRANSLATIONALLY_REGULATED_GENE = "post_translationally_regulated_gene"
RECOMBINATIONALLY_REARRANGED_VERTEBRATE_IMMUNE_SYSTEM_GENE = "recombinationally_rearranged_vertebrate_immune_system_gene"
SNORNA_GENE = "snoRNA_gene"
GENE_SILENCED_BY_HISTONE_DEACETYLATION = "gene_silenced_by_histone_deacetylation"
PROPLASTID_GENE = "proplastid_gene"
PROTEIN_CODING_GENE = "protein_coding_gene"
GENE_WITH_EDITED_TRANSCRIPT = "gene_with_edited_transcript"
RECOMBINATIONALLY_REARRANGED_GENE = "recombinationally_rearranged_gene"
LINCRNA_GENE = "lincRNA_gene"
FLOXED_GENE = "floxed_gene"
ENDOGENOUS_RETROVIRAL_GENE = "endogenous_retroviral_gene"
GENE_WITH_NON_CANONICAL_START_CODON = "gene_with_non_canonical_start_codon"
MIRNA_GENE = "miRNA_gene"
NCRNA_GENE = "ncRNA_gene"
GENE_WITH_DICISTRONIC_TRANSCRIPT = "gene_with_dicistronic_transcript"
TRANSPOSABLE_ELEMENT_GENE = "transposable_element_gene"
MT_GENE = "mt_gene"
ENGINEERED_FOREIGN_GENE = "engineered_foreign_gene"
POSITIVELY_AUTOREGULATED_GENE = "positively_autoregulated_gene"
GENE_WITH_MRNA_RECODED_BY_TRANSLATIONAL_BYPASS = "gene_with_mRNA_recoded_by_translational_bypass"
GENE_SILENCED_BY_DNA_METHYLATION = "gene_silenced_by_DNA_methylation"
LEUCOPLAST_GENE = "leucoplast_gene"
CYANELLE_GENE = "cyanelle_gene"
TMRNA_GENE = "tmRNA_gene"
MAXICIRCLE_GENE = "maxicircle_gene"

SOSuperLocusFeatures = make_enum("SOSuperLocusFeatures",
                                 GENE,
                                 EPIGENETICALLY_MODIFIED_GENE,
                                 GENE_CASSETTE,
                                 GENE_WITH_START_CODON_CUG,
                                 GENE_WITH_TRANS_SPLICED_TRANSCRIPT,
                                 TRANSLATIONALLY_REGULATED_GENE,
                                 CT_GENE,
                                 GENE_WITH_STOP_CODON_REDEFINED_AS_SELENOCYSTEINE,
                                 RECOMBINATIONALLY_INVERTED_GENE,
                                 PLASTID_GENE,
                                 GENE_WITH_DICISTRONIC_MRNA,
                                 FOREIGN_GENE,
                                 GENE_WITH_RECODED_MRNA,
                                 CRYPTOGENE,
                                 KINETOPLAST_GENE,
                                 APICOPLAST_GENE,
                                 CHROMOPLAST_GENE,
                                 GENE_WITH_STOP_CODON_READ_THROUGH,
                                 GENE_WITH_TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT,
                                 GENE_REARRANGED_AT_DNA_LEVEL,
                                 NUCLEOMORPH_GENE,
                                 TRANSGENE,
                                 FUNCTIONAL_CANDIDATE_GENE,
                                 GENE_WITH_DICISTRONIC_PRIMARY_TRANSCRIPT,
                                 PIRNA_GENE,
                                 ENGINEERED_GENE,
                                 POSITIONAL_CANDIDATE_GENE,
                                 CANDIDATE_GENE,
                                 RNASE_MRP_RNA_GENE,
                                 GENE_WITH_STOP_CODON_REDEFINED_AS_PYRROLYSINE,
                                 PROVIRAL_GENE,
                                 SRP_RNA_GENE,
                                 TRNA_GENE,
                                 CRYPTIC_GENE,
                                 FUSION_GENE,
                                 ALLELICALLY_EXCLUDED_GENE,
                                 RETROGENE,
                                 GRNA_GENE,
                                 ENGINEERED_FUSION_GENE,
                                 GENE_WITH_POLYCISTRONIC_TRANSCRIPT,
                                 PLASMID_GENE,
                                 MATERNALLY_IMPRINTED_GENE,
                                 GENE_SILENCED_BY_RNA_INTERFERENCE,
                                 WILD_TYPE_RESCUE_GENE,
                                 RNASE_P_RNA_GENE,
                                 RRNA_GENE,
                                 GENE_SILENCED_BY_DNA_MODIFICATION,
                                 ENGINEERED_FOREIGN_TRANSPOSABLE_ELEMENT_GENE,
                                 SILENCED_GENE,
                                 GENE_WITH_POLYADENYLATED_MRNA,
                                 MINICIRCLE_GENE,
                                 PREDICTED_GENE,
                                 TELOMERASE_RNA_GENE,
                                 RESCUE_GENE,
                                 PATERNALLY_IMPRINTED_GENE,
                                 NUCLEAR_GENE,
                                 GENE_WITH_MRNA_WITH_FRAMESHIFT,
                                 SCRNA_GENE,
                                 NEGATIVELY_AUTOREGULATED_GENE,
                                 GENE_SILENCED_BY_HISTONE_METHYLATION,
                                 GENE_SILENCED_BY_HISTONE_MODIFICATION,
                                 SNRNA_GENE,
                                 POST_TRANSLATIONALLY_REGULATED_GENE,
                                 RECOMBINATIONALLY_REARRANGED_VERTEBRATE_IMMUNE_SYSTEM_GENE,
                                 SNORNA_GENE,
                                 GENE_SILENCED_BY_HISTONE_DEACETYLATION,
                                 PROPLASTID_GENE,
                                 PROTEIN_CODING_GENE,
                                 GENE_WITH_EDITED_TRANSCRIPT,
                                 RECOMBINATIONALLY_REARRANGED_GENE,
                                 LINCRNA_GENE,
                                 FLOXED_GENE,
                                 ENDOGENOUS_RETROVIRAL_GENE,
                                 GENE_WITH_NON_CANONICAL_START_CODON,
                                 MIRNA_GENE,
                                 NCRNA_GENE,
                                 GENE_WITH_DICISTRONIC_TRANSCRIPT,
                                 TRANSPOSABLE_ELEMENT_GENE,
                                 MT_GENE,
                                 ENGINEERED_FOREIGN_GENE,
                                 POSITIVELY_AUTOREGULATED_GENE,
                                 GENE_WITH_MRNA_RECODED_BY_TRANSLATIONAL_BYPASS,
                                 GENE_SILENCED_BY_DNA_METHYLATION,
                                 LEUCOPLAST_GENE,
                                 CYANELLE_GENE,
                                 TMRNA_GENE,
                                 MAXICIRCLE_GENE,
                                 PSEUDOGENE,
                                 C_D_H_ACA_BOX_SCARNA_GENE,
                                 C_D_BOX_SCARNA_GENE,
                                 C_D_BOX_SNORNA_GENE,
                                 H_ACA_BOX_SCARNA_GENE,
                                 H_ACA_BOX_SNORNA_GENE,
                                 IG_C_GENE,
                                 IG_D_GENE,
                                 IG_J_GENE,
                                 IG_V_GENE,
                                 RNA_7SK_GENE,
                                 TR_C_GENE,
                                 TR_D_GENE,
                                 TR_J_GENE,
                                 TR_V_GENE,
                                 T_CELL_RECEPTOR_GENE,
                                 U14_SNORNA_GENE,
                                 U3_SNORNA_GENE,
                                 Y_RNA_GENE,
                                 ANTISENSE_LNCRNA_GENE,
                                 BIDIRECTIONAL_PROMOTER_LNCRNA_GENE,
                                 CYTOSOLIC_LSU_RRNA_GENE,
                                 CYTOSOLIC_SSU_RRNA_GENE,
                                 CYTOSOLIC_RRNA_16S_GENE,
                                 CYTOSOLIC_RRNA_18S_GENE,
                                 CYTOSOLIC_RRNA_23S_GENE,
                                 CYTOSOLIC_RRNA_25S_GENE,
                                 CYTOSOLIC_RRNA_28S_GENE,
                                 CYTOSOLIC_RRNA_2S_GENE,
                                 CYTOSOLIC_RRNA_5S_GENE,
                                 CYTOSOLIC_RRNA_5_8S_GENE,
                                 CYTOSOLIC_RRNA_GENE,
                                 ENZYMATIC_RNA_GENE,
                                 HPRNA_GENE,
                                 IMMUNOGLOBULIN_GENE,
                                 LNCRNA_GENE,
                                 METHYLATION_GUIDE_SNORNA_GENE,
                                 MT_LSU_RRNA_GENE,
                                 MT_SSU_RRNA_GENE,
                                 MT_RRNA_GENE,
                                 PLASTID_LSU_RRNA_GENE,
                                 PLASTID_SSU_RRNA_GENE,
                                 PLASTID_RRNA_GENE,
                                 PSEUDOURIDYLATION_GUIDE_SNORNA_GENE,
                                 RIBOZYME_GENE,
                                 SBRNA_GENE,
                                 SCARNA_GENE,
                                 SENSE_INTRONIC_LNCRNA_GENE,
                                 SENSE_OVERLAP_LNCRNA_GENE,
                                 SNCRNA_GENE,
                                 VAULT_RNA_GENE,
                                 VERTEBRATE_IMMUNE_SYSTEM_GENE,
)

SRP_RNA = "SRP_RNA"
HISTIDYL_TRNA = "histidyl_tRNA"
LYSYL_TRNA = "lysyl_tRNA"
U2_SNRNA = "U2_snRNA"
SERYL_TRNA = "seryl_tRNA"
ASPARTYL_TRNA = "aspartyl_tRNA"
POLYCISTRONIC_TRANSCRIPT = "polycistronic_transcript"
MICF_RNA = "MicF_RNA"
C_D_BOX_SNORNA_PRIMARY_TRANSCRIPT = "C_D_box_snoRNA_primary_transcript"
METHYLATION_GUIDE_SNORNA = "methylation_guide_snoRNA"
POLYCISTRONIC_MRNA = "polycistronic_mRNA"
MRNA_RECODED_BY_CODON_REDEFINITION = "mRNA_recoded_by_codon_redefinition"
RRNA_23S = "rRNA_23S"
U1_SNRNA = "U1_snRNA"
RASIRNA = "rasiRNA"
TYROSINE_TRNA_PRIMARY_TRANSCRIPT = "tyrosine_tRNA_primary_transcript"
TRANSCRIPT_BOUND_BY_NUCLEIC_ACID = "transcript_bound_by_nucleic_acid"
SHRNA_PRIMARY_TRANSCRIPT = "shRNA_primary_transcript"
GLUTAMINYL_TRNA = "glutaminyl_tRNA"
GLYCINE_TRNA_PRIMARY_TRANSCRIPT = "glycine_tRNA_primary_transcript"
PSEUDOURIDYLATION_GUIDE_SNORNA = "pseudouridylation_guide_snoRNA"
RRNA_25S = "rRNA_25S"
MONOCISTRONIC_PRIMARY_TRANSCRIPT = "monocistronic_primary_transcript"
U3_SNORNA = "U3_snoRNA"
LINCRNA = "lincRNA"
ANTI_ARRET = "anti_ARRET"
RPRA_RNA = "RprA_RNA"
ISOLEUCINE_TRNA_PRIMARY_TRANSCRIPT = "isoleucine_tRNA_primary_transcript"
SIRNA = "siRNA"
Y_RNA_PRIMARY_TRANSCRIPT = "Y_RNA_primary_transcript"
Y_RNA = "Y_RNA"
MIR_ENCODING_Y_RNA_PRIMARY_TRANSCRIPT = "miR_encoding_Y_RNA_primary_transcript"
PHENYLALANYL_TRNA = "phenylalanyl_tRNA"
MRNA_WITH_PLUS_2_FRAMESHIFT = "mRNA_with_plus_2_frameshift"
OXYS_RNA = "OxyS_RNA"
GUIDE_RNA = "guide_RNA"
EDITED_TRANSCRIPT_BY_A_TO_I_SUBSTITUTION = "edited_transcript_by_A_to_I_substitution"
ALTERNATIVELY_SPLICED_TRANSCRIPT = "alternatively_spliced_transcript"
TRANS_SPLICED_TRANSCRIPT = "trans_spliced_transcript"
MIR_ENCODING_LNCRNA_PRIMARY_TRANSCRIPT = "miR_encoding_lncRNA_primary_transcript"
SNORNA = "snoRNA"
PRIRNA = "priRNA"
CLASS_II_RNA = "class_II_RNA"
CLASS_I_RNA = "class_I_RNA"
RRNA_LARGE_SUBUNIT_PRIMARY_TRANSCRIPT = "rRNA_large_subunit_primary_transcript"
RRNA_5_8S = "rRNA_5_8S"
SMALL_REGULATORY_NCRNA = "small_regulatory_ncRNA"
PIRNA = "piRNA"
LEUCYL_TRNA = "leucyl_tRNA"
EDITED_TRANSCRIPT = "edited_transcript"
PYRROLYSYL_TRNA = "pyrrolysyl_tRNA"
DICISTRONIC_TRANSCRIPT = "dicistronic_transcript"
SELENOCYSTEINYL_TRNA = "selenocysteinyl_tRNA"
NCRNA = "ncRNA"
MIRNA = "miRNA"
U5_SNRNA = "U5_snRNA"
CAPPED_PRIMARY_TRANSCRIPT = "capped_primary_transcript"
C_D_BOX_SNORNA = "C_D_box_snoRNA"
PYRROLYSINE_TRNA_PRIMARY_TRANSCRIPT = "pyrrolysine_tRNA_primary_transcript"
RRNA_CLEAVAGE_SNORNA_PRIMARY_TRANSCRIPT = "rRNA_cleavage_snoRNA_primary_transcript"
LOR = "loR"
ABERRANT_PROCESSED_TRANSCRIPT = "aberrant_processed_transcript"
PROLINE_TRNA_PRIMARY_TRANSCRIPT = "proline_tRNA_primary_transcript"
MONOCISTRONIC_TRANSCRIPT = "monocistronic_transcript"
INTRONIC_LNCRNA = "intronic_lncRNA"
TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT = "transcript_with_translational_frameshift"
LARGE_SUBUNIT_RRNA = "large_subunit_rRNA"
METHYLATION_GUIDE_SNORNA_PRIMARY_TRANSCRIPT = "methylation_guide_snoRNA_primary_transcript"
GLUTAMYL_TRNA = "glutamyl_tRNA"
REGIONAL_CENTROMERE_OUTER_REPEAT_TRANSCRIPT = "regional_centromere_outer_repeat_transcript"
RRNA_28S = "rRNA_28S"
TASIRNA_PRIMARY_TRANSCRIPT = "tasiRNA_primary_transcript"
CSRB_RSMB_RNA = "CsrB_RsmB_RNA"
RRNA_SMALL_SUBUNIT_PRIMARY_TRANSCRIPT = "rRNA_small_subunit_primary_transcript"
DICISTRONIC_MRNA = "dicistronic_mRNA"
TMRNA = "tmRNA"
ARGININE_TRNA_PRIMARY_TRANSCRIPT = "arginine_tRNA_primary_transcript"
PHENYLALANINE_TRNA_PRIMARY_TRANSCRIPT = "phenylalanine_tRNA_primary_transcript"
MINI_EXON_DONOR_RNA = "mini_exon_donor_RNA"
MATURE_TRANSCRIPT = "mature_transcript"
TELOMERASE_RNA = "telomerase_RNA"
ASPARTIC_ACID_TRNA_PRIMARY_TRANSCRIPT = "aspartic_acid_tRNA_primary_transcript"
MRNA_WITH_PLUS_1_FRAMESHIFT = "mRNA_with_plus_1_frameshift"
RRNA_5S = "rRNA_5S"
RRE_RNA = "RRE_RNA"
MRNA_WITH_FRAMESHIFT = "mRNA_with_frameshift"
ANTISENSE_PRIMARY_TRANSCRIPT = "antisense_primary_transcript"
ARGINYL_TRNA = "arginyl_tRNA"
ISOLEUCYL_TRNA = "isoleucyl_tRNA"
GLUTAMIC_ACID_TRNA_PRIMARY_TRANSCRIPT = "glutamic_acid_tRNA_primary_transcript"
MRNA = "mRNA"
TRANSCRIPT = "transcript"
HISTIDINE_TRNA_PRIMARY_TRANSCRIPT = "histidine_tRNA_primary_transcript"
SRP_RNA_PRIMARY_TRANSCRIPT = "SRP_RNA_primary_transcript"
RIBOZYME = "ribozyme"
LNCRNA_PRIMARY_TRANSCRIPT = "lncRNA_primary_transcript"
METHIONYL_TRNA = "methionyl_tRNA"
ENZYMATIC_RNA = "enzymatic_RNA"
GLUTAMINE_TRNA_PRIMARY_TRANSCRIPT = "glutamine_tRNA_primary_transcript"
SPOT_42_RNA = "spot_42_RNA"
MOR = "moR"
TMRNA_PRIMARY_TRANSCRIPT = "tmRNA_primary_transcript"
GCVB_RNA = "GcvB_RNA"
U6_SNRNA = "U6_snRNA"
ASPARAGINE_TRNA_PRIMARY_TRANSCRIPT = "asparagine_tRNA_primary_transcript"
ANTISENSE_RNA = "antisense_RNA"
H_ACA_BOX_SNORNA_PRIMARY_TRANSCRIPT = "H_ACA_box_snoRNA_primary_transcript"
RRNA_21S = "rRNA_21S"
PROTEIN_CODING_PRIMARY_TRANSCRIPT = "protein_coding_primary_transcript"
SNRNA_PRIMARY_TRANSCRIPT = "snRNA_primary_transcript"
TRYPTOPHAN_TRNA_PRIMARY_TRANSCRIPT = "tryptophan_tRNA_primary_transcript"
PRIMARY_TRANSCRIPT = "primary_transcript"
NC_PRIMARY_TRANSCRIPT = "nc_primary_transcript"
U14_SNORNA_PRIMARY_TRANSCRIPT = "U14_snoRNA_primary_transcript"
U4ATAC_SNRNA = "U4atac_snRNA"
VAULTRNA_PRIMARY_TRANSCRIPT = "vaultRNA_primary_transcript"
U6ATAC_SNRNA = "U6atac_snRNA"
RNASE_P_RNA = "RNase_P_RNA"
MIR_ENCODING_TRNA_PRIMARY_TRANSCRIPT = "miR_encoding_tRNA_primary_transcript"
VALYL_TRNA = "valyl_tRNA"
EXEMPLAR_MRNA = "exemplar_mRNA"
VAULT_RNA = "vault_RNA"
SERINE_TRNA_PRIMARY_TRANSCRIPT = "serine_tRNA_primary_transcript"
POLYCISTRONIC_PRIMARY_TRANSCRIPT = "polycistronic_primary_transcript"
U12_SNRNA = "U12_snRNA"
ARRET = "ARRET"
SMALL_SUBUNIT_RRNA = "small_subunit_rRNA"
U14_SNORNA = "U14_snoRNA"
THREONINE_TRNA_PRIMARY_TRANSCRIPT = "threonine_tRNA_primary_transcript"
RRNA = "rRNA"
CYSTEINE_TRNA_PRIMARY_TRANSCRIPT = "cysteine_tRNA_primary_transcript"
H_ACA_BOX_SNORNA = "H_ACA_box_snoRNA"
METHIONINE_TRNA_PRIMARY_TRANSCRIPT = "methionine_tRNA_primary_transcript"
CONSENSUS_MRNA = "consensus_mRNA"
U11_SNRNA = "U11_snRNA"
RNASE_MRP_RNA = "RNase_MRP_RNA"
MRNA_WITH_MINUS_1_FRAMESHIFT = "mRNA_with_minus_1_frameshift"
TELOMERIC_TRANSCRIPT = "telomeric_transcript"
ENHANCERRNA = "enhancerRNA"
RNA_6S = "RNA_6S"
TERRA = "TERRA"
LYSINE_TRNA_PRIMARY_TRANSCRIPT = "lysine_tRNA_primary_transcript"
CYSTEINYL_TRNA = "cysteinyl_tRNA"
SNORNA_PRIMARY_TRANSCRIPT = "snoRNA_primary_transcript"
SNRNA = "snRNA"
TRANS_SPLICED_MRNA = "trans_spliced_mRNA"
DSRA_RNA = "DsrA_RNA"
RRNA_18S = "rRNA_18S"
MIR_ENCODING_SNORNA_PRIMARY_TRANSCRIPT = "miR_encoding_snoRNA_primary_transcript"
MIR_ENCODING_SHRNA_PRIMARY_TRANSCRIPT = "miR_encoding_shRNA_primary_transcript"
CAPPED_MRNA = "capped_mRNA"
SELENOCYSTEINE_TRNA_PRIMARY_TRANSCRIPT = "selenocysteine_tRNA_primary_transcript"
EDITED_MRNA = "edited_mRNA"
ALANINE_TRNA_PRIMARY_TRANSCRIPT = "alanine_tRNA_primary_transcript"
SHRNA = "shRNA"
RECODED_MRNA = "recoded_mRNA"
TRNA = "tRNA"
TRNA_PRIMARY_TRANSCRIPT = "tRNA_primary_transcript"
ASPARAGINYL_TRNA = "asparaginyl_tRNA"
THREONYL_TRNA = "threonyl_tRNA"
TYROSYL_TRNA = "tyrosyl_tRNA"
ANTISENSE_LNCRNA = "antisense_lncRNA"
POLYADENYLATED_MRNA = "polyadenylated_mRNA"
PROLYL_TRNA = "prolyl_tRNA"
ARIA = "ARIA"
SCRNA = "scRNA"
PROCESSED_TRANSCRIPT = "processed_transcript"
U4_SNRNA = "U4_snRNA"
LEUCINE_TRNA_PRIMARY_TRANSCRIPT = "leucine_tRNA_primary_transcript"
ALANYL_TRNA = "alanyl_tRNA"
GLYCYL_TRNA = "glycyl_tRNA"
VALINE_TRNA_PRIMARY_TRANSCRIPT = "valine_tRNA_primary_transcript"
MRNA_WITH_MINUS_2_FRAMESHIFT = "mRNA_with_minus_2_frameshift"
MIR_ENCODING_VAULTRNA_PRIMARY_TRANSCRIPT = "miR_encoding_vaultRNA_primary_transcript"
TRYPTOPHANYL_TRNA = "tryptophanyl_tRNA"
SCRNA_PRIMARY_TRANSCRIPT = "scRNA_primary_transcript"
RRNA_CLEAVAGE_RNA = "rRNA_cleavage_RNA"
MRNA_RECODED_BY_TRANSLATIONAL_BYPASS = "mRNA_recoded_by_translational_bypass"
RRNA_16S = "rRNA_16S"
TRANSCRIPT_BOUND_BY_PROTEIN = "transcript_bound_by_protein"
PRE_EDITED_MRNA = "pre_edited_mRNA"
DICISTRONIC_PRIMARY_TRANSCRIPT = "dicistronic_primary_transcript"
MIRNA_PRIMARY_TRANSCRIPT = "miRNA_primary_transcript"
MONOCISTRONIC_MRNA = "monocistronic_mRNA"
RRNA_PRIMARY_TRANSCRIPT = "rRNA_primary_transcript"
TASIRNA = "tasiRNA"

SOTranscriptFeatures = make_enum("SOTranscriptFeatures",
                                 SRP_RNA,
                                 HISTIDYL_TRNA,
                                 LYSYL_TRNA,
                                 U2_SNRNA,
                                 SERYL_TRNA,
                                 ASPARTYL_TRNA,
                                 POLYCISTRONIC_TRANSCRIPT,
                                 MICF_RNA,
                                 C_D_BOX_SNORNA_PRIMARY_TRANSCRIPT,
                                 METHYLATION_GUIDE_SNORNA,
                                 POLYCISTRONIC_MRNA,
                                 MRNA_RECODED_BY_CODON_REDEFINITION,
                                 RRNA_23S,
                                 U1_SNRNA,
                                 RASIRNA,
                                 TYROSINE_TRNA_PRIMARY_TRANSCRIPT,
                                 TRANSCRIPT_BOUND_BY_NUCLEIC_ACID,
                                 SHRNA_PRIMARY_TRANSCRIPT,
                                 GLUTAMINYL_TRNA,
                                 GLYCINE_TRNA_PRIMARY_TRANSCRIPT,
                                 PSEUDOURIDYLATION_GUIDE_SNORNA,
                                 RRNA_25S,
                                 MONOCISTRONIC_PRIMARY_TRANSCRIPT,
                                 U3_SNORNA,
                                 LINCRNA,
                                 ANTI_ARRET,
                                 RPRA_RNA,
                                 ISOLEUCINE_TRNA_PRIMARY_TRANSCRIPT,
                                 SIRNA,
                                 Y_RNA_PRIMARY_TRANSCRIPT,
                                 Y_RNA,
                                 MIR_ENCODING_Y_RNA_PRIMARY_TRANSCRIPT,
                                 PHENYLALANYL_TRNA,
                                 MRNA_WITH_PLUS_2_FRAMESHIFT,
                                 OXYS_RNA,
                                 GUIDE_RNA,
                                 EDITED_TRANSCRIPT_BY_A_TO_I_SUBSTITUTION,
                                 ALTERNATIVELY_SPLICED_TRANSCRIPT,
                                 TRANS_SPLICED_TRANSCRIPT,
                                 MIR_ENCODING_LNCRNA_PRIMARY_TRANSCRIPT,
                                 SNORNA,
                                 PRIRNA,
                                 CLASS_II_RNA,
                                 CLASS_I_RNA,
                                 RRNA_LARGE_SUBUNIT_PRIMARY_TRANSCRIPT,
                                 RRNA_5_8S,
                                 SMALL_REGULATORY_NCRNA,
                                 PIRNA,
                                 LEUCYL_TRNA,
                                 EDITED_TRANSCRIPT,
                                 PYRROLYSYL_TRNA,
                                 DICISTRONIC_TRANSCRIPT,
                                 SELENOCYSTEINYL_TRNA,
                                 NCRNA,
                                 MIRNA,
                                 U5_SNRNA,
                                 CAPPED_PRIMARY_TRANSCRIPT,
                                 C_D_BOX_SNORNA,
                                 PYRROLYSINE_TRNA_PRIMARY_TRANSCRIPT,
                                 RRNA_CLEAVAGE_SNORNA_PRIMARY_TRANSCRIPT,
                                 LOR,
                                 ABERRANT_PROCESSED_TRANSCRIPT,
                                 PROLINE_TRNA_PRIMARY_TRANSCRIPT,
                                 MONOCISTRONIC_TRANSCRIPT,
                                 INTRONIC_LNCRNA,
                                 TRANSCRIPT_WITH_TRANSLATIONAL_FRAMESHIFT,
                                 LARGE_SUBUNIT_RRNA,
                                 METHYLATION_GUIDE_SNORNA_PRIMARY_TRANSCRIPT,
                                 GLUTAMYL_TRNA,
                                 REGIONAL_CENTROMERE_OUTER_REPEAT_TRANSCRIPT,
                                 RRNA_28S,
                                 TASIRNA_PRIMARY_TRANSCRIPT,
                                 CSRB_RSMB_RNA,
                                 RRNA_SMALL_SUBUNIT_PRIMARY_TRANSCRIPT,
                                 DICISTRONIC_MRNA,
                                 TMRNA,
                                 ARGININE_TRNA_PRIMARY_TRANSCRIPT,
                                 PHENYLALANINE_TRNA_PRIMARY_TRANSCRIPT,
                                 MINI_EXON_DONOR_RNA,
                                 MATURE_TRANSCRIPT,
                                 TELOMERASE_RNA,
                                 ASPARTIC_ACID_TRNA_PRIMARY_TRANSCRIPT,
                                 MRNA_WITH_PLUS_1_FRAMESHIFT,
                                 RRNA_5S,
                                 RRE_RNA,
                                 MRNA_WITH_FRAMESHIFT,
                                 ANTISENSE_PRIMARY_TRANSCRIPT,
                                 ARGINYL_TRNA,
                                 ISOLEUCYL_TRNA,
                                 GLUTAMIC_ACID_TRNA_PRIMARY_TRANSCRIPT,
                                 MRNA,
                                 TRANSCRIPT,
                                 HISTIDINE_TRNA_PRIMARY_TRANSCRIPT,
                                 SRP_RNA_PRIMARY_TRANSCRIPT,
                                 RIBOZYME,
                                 LNCRNA_PRIMARY_TRANSCRIPT,
                                 METHIONYL_TRNA,
                                 ENZYMATIC_RNA,
                                 GLUTAMINE_TRNA_PRIMARY_TRANSCRIPT,
                                 SPOT_42_RNA,
                                 MOR,
                                 TMRNA_PRIMARY_TRANSCRIPT,
                                 GCVB_RNA,
                                 U6_SNRNA,
                                 ASPARAGINE_TRNA_PRIMARY_TRANSCRIPT,
                                 ANTISENSE_RNA,
                                 H_ACA_BOX_SNORNA_PRIMARY_TRANSCRIPT,
                                 RRNA_21S,
                                 PROTEIN_CODING_PRIMARY_TRANSCRIPT,
                                 SNRNA_PRIMARY_TRANSCRIPT,
                                 TRYPTOPHAN_TRNA_PRIMARY_TRANSCRIPT,
                                 PRIMARY_TRANSCRIPT,
                                 NC_PRIMARY_TRANSCRIPT,
                                 U14_SNORNA_PRIMARY_TRANSCRIPT,
                                 U4ATAC_SNRNA,
                                 VAULTRNA_PRIMARY_TRANSCRIPT,
                                 U6ATAC_SNRNA,
                                 LNC_RNA,
                                 RNASE_P_RNA,
                                 MIR_ENCODING_TRNA_PRIMARY_TRANSCRIPT,
                                 VALYL_TRNA,
                                 EXEMPLAR_MRNA,
                                 VAULT_RNA,
                                 SERINE_TRNA_PRIMARY_TRANSCRIPT,
                                 POLYCISTRONIC_PRIMARY_TRANSCRIPT,
                                 U12_SNRNA,
                                 ARRET,
                                 SMALL_SUBUNIT_RRNA,
                                 U14_SNORNA,
                                 THREONINE_TRNA_PRIMARY_TRANSCRIPT,
                                 RRNA,
                                 CYSTEINE_TRNA_PRIMARY_TRANSCRIPT,
                                 H_ACA_BOX_SNORNA,
                                 METHIONINE_TRNA_PRIMARY_TRANSCRIPT,
                                 CONSENSUS_MRNA,
                                 U11_SNRNA,
                                 RNASE_MRP_RNA,
                                 MRNA_WITH_MINUS_1_FRAMESHIFT,
                                 TELOMERIC_TRANSCRIPT,
                                 ENHANCERRNA,
                                 RNA_6S,
                                 TERRA,
                                 LYSINE_TRNA_PRIMARY_TRANSCRIPT,
                                 CYSTEINYL_TRNA,
                                 SNORNA_PRIMARY_TRANSCRIPT,
                                 SNRNA,
                                 TRANS_SPLICED_MRNA,
                                 DSRA_RNA,
                                 RRNA_18S,
                                 MIR_ENCODING_SNORNA_PRIMARY_TRANSCRIPT,
                                 MIR_ENCODING_SHRNA_PRIMARY_TRANSCRIPT,
                                 CAPPED_MRNA,
                                 SELENOCYSTEINE_TRNA_PRIMARY_TRANSCRIPT,
                                 EDITED_MRNA,
                                 ALANINE_TRNA_PRIMARY_TRANSCRIPT,
                                 SHRNA,
                                 RECODED_MRNA,
                                 TRNA,
                                 TRNA_PRIMARY_TRANSCRIPT,
                                 ASPARAGINYL_TRNA,
                                 THREONYL_TRNA,
                                 TYROSYL_TRNA,
                                 ANTISENSE_LNCRNA,
                                 POLYADENYLATED_MRNA,
                                 PROLYL_TRNA,
                                 ARIA,
                                 SCRNA,
                                 PROCESSED_TRANSCRIPT,
                                 U4_SNRNA,
                                 LEUCINE_TRNA_PRIMARY_TRANSCRIPT,
                                 ALANYL_TRNA,
                                 GLYCYL_TRNA,
                                 VALINE_TRNA_PRIMARY_TRANSCRIPT,
                                 MRNA_WITH_MINUS_2_FRAMESHIFT,
                                 MIR_ENCODING_VAULTRNA_PRIMARY_TRANSCRIPT,
                                 TRYPTOPHANYL_TRNA,
                                 SCRNA_PRIMARY_TRANSCRIPT,
                                 RRNA_CLEAVAGE_RNA,
                                 MRNA_RECODED_BY_TRANSLATIONAL_BYPASS,
                                 RRNA_16S,
                                 TRANSCRIPT_BOUND_BY_PROTEIN,
                                 PRE_EDITED_MRNA,
                                 DICISTRONIC_PRIMARY_TRANSCRIPT,
                                 MIRNA_PRIMARY_TRANSCRIPT,
                                 MONOCISTRONIC_MRNA,
                                 RRNA_PRIMARY_TRANSCRIPT,
                                 TASIRNA,
                                 C_D_H_ACA_BOX_SCARNA,
                                 C_D_BOX_SCARNA,
                                 H_ACA_BOX_SCARNA,
                                 NMD_POLYMORPHIC_PSEUDOGENE_TRANSCRIPT,
                                 NMD_TRANSCRIPT,
                                 NSD_TRANSCRIPT,
                                 RNA_7SK,
                                 U7_SNRNA,
                                 BIDIRECTIONAL_PROMOTER_LNCRNA,
                                 CIRCULAR_MRNA,
                                 CIRCULAR_NCRNA,
                                 CODING_TRANSCRIPT_WITH_RETAINED_INTRON,
                                 CYTOSOLIC_16S_RRNA,
                                 CYTOSOLIC_18S_RRNA,
                                 CYTOSOLIC_23S_RRNA,
                                 CYTOSOLIC_25S_RRNA,
                                 CYTOSOLIC_28S_RRNA,
                                 CYTOSOLIC_2S_RRNA,
                                 CYTOSOLIC_5S_RRNA,
                                 CYTOSOLIC_5_8S_RRNA,
                                 CYTOSOLIC_LSU_RRNA,
                                 CYTOSOLIC_SSU_RRNA,
                                 CYTOSOLIC_RRNA,
                                 HPRNA,
                                 ISOMIR,
                                 LNCRNA_WITH_RETAINED_INTRON,
                                 LONG_TERMINAL_REPEAT_TRANSCRIPT,
                                 MT_LSU_RRNA,
                                 MT_SSU_RRNA,
                                 MT_RRNA,
                                 MT_TRNA,
                                 PLASTID_LSU_RRNA,
                                 PLASTID_SSU_RRNA,
                                 PLASTID_RRNA,
                                 POLYMORPHIC_PSEUDOGENE_PROCESSED_TRANSCRIPT,
                                 PRECURSOR_RNA,
                                 PREDICTED_TRANSCRIPT,
                                 PSEUDOGENE_PROCESSED_TRANSCRIPT,
                                 PSEUDOGENIC_TRANSCRIPT_WITH_RETAINED_INTRON,
                                 REF_MIRNA,
                                 SBRNA,
                                 SCARNA,
                                 SELF_CLEAVING_RIBOZYME,
                                 SENSE_INTRONIC_LNCRNA,
                                 SENSE_OVERLAP_LNCRNA,
                                 SISRNA,
                                 SNCRNA,
                                 THREE_PRIME_OVERLAPPING_NCRNA,
                                 UNCONFIRMED_TRANSCRIPT,
)

THREE_PRIME_CODING_EXON = "three_prime_coding_exon"
INTERIOR_EXON = "interior_exon"
INTERIOR_CODING_EXON = "interior_coding_exon"
EXON = "exon"
THREE_PRIME_NONCODING_EXON = "three_prime_noncoding_exon"
CODING_EXON = "coding_exon"
EXON_OF_SINGLE_EXON_GENE = "exon_of_single_exon_gene"
FIVE_PRIME_NONCODING_EXON = "five_prime_noncoding_exon"
FIVE_PRIME_CODING_EXON = "five_prime_coding_exon"
NONCODING_EXON = "noncoding_exon"

SOExonFeatures = make_enum("SOExonFeatures",
                           THREE_PRIME_CODING_EXON,
                           INTERIOR_EXON,
                           INTERIOR_CODING_EXON,
                           EXON,
                           THREE_PRIME_NONCODING_EXON,
                           CODING_EXON,
                           EXON_OF_SINGLE_EXON_GENE,
                           FIVE_PRIME_NONCODING_EXON,
                           FIVE_PRIME_CODING_EXON,
                           NONCODING_EXON,
)

CDS_SUPPORTED_BY_EST_OR_CDNA_DATA = "CDS_supported_by_EST_or_cDNA_data"
CDS_PREDICTED = "CDS_predicted"
CDS_SUPPORTED_BY_DOMAIN_MATCH_DATA = "CDS_supported_by_domain_match_data"
CDS_INDEPENDENTLY_KNOWN = "CDS_independently_known"
CDS_SUPPORTED_BY_SEQUENCE_SIMILARITY_DATA = "CDS_supported_by_sequence_similarity_data"
EDITED_CDS = "edited_CDS"
TRANSPOSABLE_ELEMENT_CDS = "transposable_element_CDS"
ORPHAN_CDS = "orphan_CDS"
CDS_FRAGMENT = "CDS_fragment"
CDS = "CDS"

SOCDSFeatures = make_enum("SOCDSFeatures",
                          CDS_SUPPORTED_BY_EST_OR_CDNA_DATA,
                          CDS_PREDICTED,
                          CDS_SUPPORTED_BY_DOMAIN_MATCH_DATA,
                          CDS_INDEPENDENTLY_KNOWN,
                          CDS_SUPPORTED_BY_SEQUENCE_SIMILARITY_DATA,
                          EDITED_CDS,
                          TRANSPOSABLE_ELEMENT_CDS,
                          ORPHAN_CDS,
                          CDS_FRAGMENT,
                          CDS,
                          CDS_EXTENSION,
                          CDS_FIVE_PRIME_EXTENSION,
                          CDS_SUPPORTED_BY_PEPTIDE_SPECTRUM_MATCH,
                          CDS_THREE_PRIME_EXTENSION,
)

# the untranslated part of an exon; whole untranslated exons (e.g. noncoding_exon) are exon
# features, and UTR_region is a region within a UTR rather than a piece of one
SOUTRFeatures = make_enum("SOUTRFeatures",
                          UTR,
                          FIVE_PRIME_UTR,
                          THREE_PRIME_UTR,
                          INTERNAL_UTR,
                          NONCODING_REGION_OF_EXON,
                          FIVE_PRIME_CODING_EXON_NONCODING_REGION,
                          THREE_PRIME_CODING_EXON_NONCODING_REGION,
                          UNTRANSLATED_REGION_POLYCISTRONIC_MRNA,
)


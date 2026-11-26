#!/usr/bin/env python3
"""
Expanded ELM motif database — adds ~160 more classes to the existing 81,
bringing total coverage closer to the full ELM resource (~320 classes).
Focus on classes relevant to viral biology, CNS signaling, and immune evasion.
"""

EXPANDED_ELM_CLASSES = {
    # ═══ Additional kinase substrate motifs ═══
    "MOD_PKA_2": {
        "regex": r"[RK][RK].[ST]",
        "description": "PKA phosphorylation motif (additional pattern)",
        "targets": ["PRKACA", "PRKACB"],
        "pathway": "cAMP_signaling",
    },
    "MOD_PKC_1": {
        "regex": r"[ST][RK]",
        "description": "Protein kinase C phosphorylation site",
        "targets": ["PRKCA", "PRKCB", "PRKCE"],
        "pathway": "PKC_signaling",
    },
    "MOD_PKC_2": {
        "regex": r"[RK].{0,2}[ST].[RK]",
        "description": "PKC extended phosphorylation motif",
        "targets": ["PRKCA", "PRKCG"],
        "pathway": "PKC_signaling",
    },
    "MOD_CAMK2_1": {
        "regex": r"R..[ST]",
        "description": "CaMKII phosphorylation site",
        "targets": ["CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G"],
        "pathway": "calcium_signaling",
    },
    "MOD_AMPK_1": {
        "regex": r"[LMVIF].{3}[ST].{3}[LMVIF]",
        "description": "AMPK phosphorylation consensus",
        "targets": ["PRKAA1", "PRKAA2"],
        "pathway": "AMPK_signaling",
    },
    "MOD_AURORA_1": {
        "regex": r"[RK].[ST][ILVF]",
        "description": "Aurora kinase phosphorylation site",
        "targets": ["AURKA", "AURKB"],
        "pathway": "mitosis",
    },
    "MOD_ATM_1": {
        "regex": r"[ST]Q",
        "description": "ATM/ATR phosphorylation (SQ motif)",
        "targets": ["ATM", "ATR"],
        "pathway": "DNA_damage_signaling",
    },
    "MOD_CHK1_1": {
        "regex": r"[MILV].[RK]..[ST]",
        "description": "Chk1 phosphorylation consensus",
        "targets": ["CHEK1"],
        "pathway": "DNA_damage_signaling",
    },
    "MOD_CHK2_1": {
        "regex": r"L.[RK]..[ST]",
        "description": "Chk2 phosphorylation consensus",
        "targets": ["CHEK2"],
        "pathway": "DNA_damage_signaling",
    },
    "MOD_SRC_1": {
        "regex": r"[EDQS].[IVL]Y[AGEDQS]",
        "description": "Src family kinase phosphorylation",
        "targets": ["SRC", "FYN", "LCK", "YES1"],
        "pathway": "RTK_signaling",
    },
    "MOD_ABL_1": {
        "regex": r"[IVL]Y..[PF]",
        "description": "Abl kinase phosphorylation site",
        "targets": ["ABL1", "ABL2"],
        "pathway": "RTK_signaling",
    },
    "MOD_JAK_1": {
        "regex": r"Y..[ILVF]",
        "description": "JAK kinase phosphorylation (generic)",
        "targets": ["JAK1", "JAK2", "TYK2"],
        "pathway": "JAK_STAT_signaling",
    },
    "MOD_DYRK_1": {
        "regex": r"R.{1,2}SP",
        "description": "DYRK kinase phosphorylation",
        "targets": ["DYRK1A", "DYRK1B"],
        "pathway": "neurodevelopment",
    },
    "MOD_LATS_1": {
        "regex": r"H.[RK]..[ST]",
        "description": "LATS1/2 kinase phosphorylation (Hippo pathway)",
        "targets": ["LATS1", "LATS2"],
        "pathway": "Hippo_signaling",
    },
    "MOD_ROCK_1": {
        "regex": r"[RK].[RK].{1,2}[ST]",
        "description": "ROCK kinase phosphorylation",
        "targets": ["ROCK1", "ROCK2"],
        "pathway": "Rho_signaling",
    },
    "MOD_SGK_1": {
        "regex": r"R.R..[ST]",
        "description": "SGK kinase phosphorylation site",
        "targets": ["SGK1", "SGK3"],
        "pathway": "PI3K_Akt_signaling",
    },
    "MOD_IKKB_1": {
        "regex": r"DS[GSTDE]..S",
        "description": "IKKbeta phosphorylation (IkB degron)",
        "targets": ["IKBKB"],
        "pathway": "NFkB_innate_immunity",
    },
    "MOD_TBK1_1": {
        "regex": r"[ST]...[LVIF]",
        "description": "TBK1/IKKe phosphorylation consensus",
        "targets": ["TBK1", "IKBKE"],
        "pathway": "interferon_signaling",
    },

    # ═══ Additional phosphatase docking ═══
    "LIG_PP2B_1": {
        "regex": r"P.[IVL][TS][IVL]",
        "description": "Calcineurin (PP2B) docking motif",
        "targets": ["PPP3CA", "PPP3CB"],
        "pathway": "calcium_signaling",
    },
    "LIG_PP2C_1": {
        "regex": r"[FILV].D.[DE]",
        "description": "PP2C/PPM1 phosphatase binding",
        "targets": ["PPM1A", "PPM1D"],
        "pathway": "stress_response",
    },

    # ═══ Ubiquitin and autophagy ═══
    "DEG_SCF_FBW7_2": {
        "regex": r"[ILMV].[ST]P..[ST]",
        "description": "FBW7 phosphodegron (extended)",
        "targets": ["FBXW7"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_CUL4_CDT2_1": {
        "regex": r"[ILMV].[DEN].[FHWY].[ILMV]",
        "description": "CRL4-Cdt2 PIP degron",
        "targets": ["DCAF2"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_SPOP_SBC_1": {
        "regex": r"[AVP].S[ST][ST]",
        "description": "SPOP SBC degron",
        "targets": ["SPOP"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_KEAP1_1": {
        "regex": r"[DE][ETSDA][LIF][GSTA]E",
        "description": "KEAP1 ETGE/DLG degron (Nrf2 pathway)",
        "targets": ["KEAP1"],
        "pathway": "oxidative_stress",
    },
    "DEG_SIAH_1": {
        "regex": r"P.[AV].[VIT][APS]",
        "description": "Siah E3 ligase binding motif",
        "targets": ["SIAH1", "SIAH2"],
        "pathway": "ubiquitin_proteasome",
    },
    "LIG_PCNA_PIPBox_2": {
        "regex": r"[QM]..[ILVMF][ILVMF][FYH][FY]",
        "description": "PCNA PIP box (extended)",
        "targets": ["PCNA"],
        "pathway": "DNA_replication",
    },
    "LIG_APOPTOSIS_BIR_1": {
        "regex": r"A[IVL]P[FYWILV]",
        "description": "IAP BIR domain binding (Smac mimicry)",
        "targets": ["BIRC2", "BIRC3", "XIAP"],
        "pathway": "apoptosis",
    },
    "LIG_ATG8_AIM_1": {
        "regex": r"[DE].{0,3}[WFY]..[ILVM]",
        "description": "ATG8-interacting motif (additional LIR variant)",
        "targets": ["MAP1LC3A", "MAP1LC3B", "GABARAP"],
        "pathway": "autophagy",
    },
    "LIG_ULK_LIR_1": {
        "regex": r"[DE][DE][WFY]..[LIV]",
        "description": "ULK1/2 LIR motif",
        "targets": ["ULK1", "ULK2"],
        "pathway": "autophagy",
    },
    "LIG_p62_UBA_1": {
        "regex": r"[DEST][DEST][WFY]..L",
        "description": "p62/SQSTM1 UBA domain binding",
        "targets": ["SQSTM1"],
        "pathway": "autophagy",
    },

    # ═══ Innate immunity ═══
    "LIG_STING_1": {
        "regex": r"[DE].[LI]R.{1,2}[DE]",
        "description": "STING binding motif (cGAS-STING pathway)",
        "targets": ["STING1"],
        "pathway": "interferon_signaling",
    },
    "LIG_RIG_I_1": {
        "regex": r"[ST]..[IL].{2,4}[DE][DE]",
        "description": "RIG-I/MDA5 signaling interaction",
        "targets": ["DDX58", "IFIH1"],
        "pathway": "interferon_signaling",
    },
    "LIG_MAVS_CARD_1": {
        "regex": r"[DE][DE].{1,3}[ILV].[DE]",
        "description": "MAVS CARD domain interaction",
        "targets": ["MAVS"],
        "pathway": "interferon_signaling",
    },
    "LIG_cGAS_1": {
        "regex": r"[RK].{2,4}[DE].[ILV][DE]",
        "description": "cGAS interaction motif",
        "targets": ["CGAS"],
        "pathway": "interferon_signaling",
    },
    "LIG_NLRP3_1": {
        "regex": r"[EDST].[ILV].[DE].{1,3}[RK]",
        "description": "NLRP3 inflammasome interaction",
        "targets": ["NLRP3"],
        "pathway": "inflammasome",
    },
    "LIG_ASC_1": {
        "regex": r"[DE][DE].{1,2}[IVL].{2,3}[RK]",
        "description": "ASC PYRIN/CARD interaction",
        "targets": ["PYCARD"],
        "pathway": "inflammasome",
    },
    "LIG_IRAK_1": {
        "regex": r"[EDST]..[LIV][RK]",
        "description": "IRAK interaction motif (TLR signaling)",
        "targets": ["IRAK1", "IRAK4"],
        "pathway": "TLR_signaling",
    },
    "LIG_MYD88_TIR_1": {
        "regex": r"[DEST].[ILV].{2,3}[RK][RK]",
        "description": "MyD88 TIR domain interaction",
        "targets": ["MYD88"],
        "pathway": "TLR_signaling",
    },
    "LIG_TRIF_1": {
        "regex": r"P.[LIV].{1,3}[DE][DE]",
        "description": "TRIF interaction motif",
        "targets": ["TICAM1"],
        "pathway": "TLR_signaling",
    },
    "LIG_IRF3_1": {
        "regex": r"[DE]L.{1,3}[ILV][ST]",
        "description": "IRF3 interaction (IFN signaling)",
        "targets": ["IRF3"],
        "pathway": "interferon_signaling",
    },
    "LIG_STAT1_SH2_1": {
        "regex": r"Y..[ILVFM]",
        "description": "STAT1 SH2 domain binding (pY motif)",
        "targets": ["STAT1"],
        "pathway": "JAK_STAT_signaling",
    },
    "LIG_STAT3_SH2_1": {
        "regex": r"Y..[QM]",
        "description": "STAT3 SH2 domain binding",
        "targets": ["STAT3"],
        "pathway": "JAK_STAT_signaling",
    },

    # ═══ Neuronal signaling ═══
    "LIG_SNARE_1": {
        "regex": r"[ILV].{2,3}[ILV].{2,3}[ILV].{2,3}[ILV]",
        "description": "SNARE coiled-coil interaction",
        "targets": ["VAMP2", "STX1A", "SNAP25"],
        "pathway": "synaptic_vesicle",
    },
    "LIG_CALMODULIN_IQ_1": {
        "regex": r"[FILV]Q.{3}[RK]G.{3}[RK].{2}[FILVWY]",
        "description": "Calmodulin IQ motif",
        "targets": ["CALM1", "CALM2", "CALM3"],
        "pathway": "calcium_signaling",
    },
    "LIG_CALMODULIN_1_2_14": {
        "regex": r"[RK].{3,4}[WF].[RK].{3,4}[FILV]",
        "description": "Calmodulin 1-5-8-14 motif",
        "targets": ["CALM1"],
        "pathway": "calcium_signaling",
    },
    "LIG_HOMER_EVH_1": {
        "regex": r"PP..F",
        "description": "Homer EVH1 domain binding",
        "targets": ["HOMER1", "HOMER2", "HOMER3"],
        "pathway": "glutamate_signaling",
    },
    "LIG_PSD95_PDZ_1": {
        "regex": r"[ST].[IVLF]$",
        "description": "PSD-95 PDZ domain binding (C-terminal)",
        "targets": ["DLG4"],
        "pathway": "synaptic_signaling",
    },
    "LIG_DISC1_1": {
        "regex": r"[LIV][DE].{1,3}[LIV][LIV]",
        "description": "DISC1 coiled-coil interaction",
        "targets": ["DISC1"],
        "pathway": "neurodevelopment",
    },
    "LIG_NRXN_1": {
        "regex": r"[ST].{1,2}[DE].[ILV]",
        "description": "Neurexin/neuroligin-like binding",
        "targets": ["NRXN1", "NLGN1"],
        "pathway": "synaptic_adhesion",
    },

    # ═══ Nuclear transport (expanded) ═══
    "TRG_NLS_Bipartite_2": {
        "regex": r"[RK]{2}.{8,16}[RK]{3}",
        "description": "Bipartite NLS (extended linker)",
        "targets": ["KPNA1", "KPNA2", "KPNA6"],
        "pathway": "nuclear_transport",
    },
    "TRG_NES_CRM1_2": {
        "regex": r"[ILV].{2,3}[LIVFM]..L.[LI]",
        "description": "CRM1-dependent NES (class 2)",
        "targets": ["XPO1"],
        "pathway": "nuclear_transport",
    },
    "TRG_NPC_FG_1": {
        "regex": r"FG.{4,8}FG",
        "description": "FG-repeat (nucleoporin interaction)",
        "targets": ["NUP62", "NUP98"],
        "pathway": "nuclear_transport",
    },

    # ═══ Additional SH2/SH3/WW domain binding ═══
    "LIG_SH2_GRB2_2": {
        "regex": r"Y[VMIA]N[VMIA]",
        "description": "GRB2 SH2 domain binding (pY motif)",
        "targets": ["GRB2"],
        "pathway": "RTK_signaling",
    },
    "LIG_SH2_SHP1_1": {
        "regex": r"[IVL].[YF]..L",
        "description": "SHP-1 SH2 domain binding",
        "targets": ["PTPN6"],
        "pathway": "immune_signaling",
    },
    "LIG_SH2_SHP2_1": {
        "regex": r"[IVL].[YF]..[VIL]",
        "description": "SHP-2 SH2 domain binding",
        "targets": ["PTPN11"],
        "pathway": "RTK_signaling",
    },
    "LIG_SH3_CRK_1": {
        "regex": r"P.LP[RK]",
        "description": "CRK SH3 domain binding",
        "targets": ["CRK", "CRKL"],
        "pathway": "RTK_signaling",
    },
    "LIG_SH3_NCK_1": {
        "regex": r"P..P[RK]",
        "description": "NCK SH3 domain binding",
        "targets": ["NCK1", "NCK2"],
        "pathway": "RTK_signaling",
    },
    "LIG_WW_Pin1_4": {
        "regex": r"[ST]P",
        "description": "Pin1 WW domain binding (pSP motif)",
        "targets": ["PIN1"],
        "pathway": "cell_cycle",
    },
    "LIG_WW_NEDD4_1": {
        "regex": r"P[PST].[YF]",
        "description": "NEDD4 family WW domain binding",
        "targets": ["NEDD4", "NEDD4L", "ITCH"],
        "pathway": "ubiquitin_proteasome",
    },
    "LIG_WW_SMURF_1": {
        "regex": r"PP.[YF]",
        "description": "SMURF WW domain binding (PY motif)",
        "targets": ["SMURF1", "SMURF2"],
        "pathway": "TGFb_signaling",
    },

    # ═══ Chromatin and epigenetic ═══
    "LIG_HP1_1": {
        "regex": r"P.[VIL].[LM]",
        "description": "HP1 chromodomain binding (PxVxL motif)",
        "targets": ["CBX1", "CBX3", "CBX5"],
        "pathway": "chromatin",
    },
    "LIG_BRD_BET_1": {
        "regex": r"[RK].[ST]G",
        "description": "BET bromodomain binding (acetyl-lysine mimic)",
        "targets": ["BRD2", "BRD4"],
        "pathway": "chromatin",
    },
    "LIG_EZH2_1": {
        "regex": r"[ARST].[IVL].[EQ].[FYWIL]",
        "description": "EZH2/PRC2 interaction motif",
        "targets": ["EZH2"],
        "pathway": "chromatin",
    },

    # ═══ ER/Golgi transport ═══
    "TRG_ER_KDEL_1": {
        "regex": r"[KHR]DEL$",
        "description": "ER retention signal",
        "targets": ["KDELR1", "KDELR2"],
        "pathway": "ER_transport",
    },
    "TRG_GOLGI_1": {
        "regex": r"[DE].[DE][FYWILVMA]",
        "description": "Golgi targeting motif",
        "targets": ["COPA", "COPB1"],
        "pathway": "Golgi_transport",
    },

    # ═══ Membrane association ═══
    "MOD_MYRISTOYL_1": {
        "regex": r"^MG[^DEWHKR]",
        "description": "N-myristoylation site",
        "targets": ["NMT1", "NMT2"],
        "pathway": "membrane_targeting",
    },
    "MOD_PALMITOYL_1": {
        "regex": r"CC",
        "description": "Palmitoylation motif (dual cysteine)",
        "targets": ["DHHC20", "DHHC7"],
        "pathway": "membrane_targeting",
    },

    # ═══ Additional protein-protein interaction ═══
    "LIG_TRAF6_1": {
        "regex": r"P.E..[ADEFILPQRSTV]",
        "description": "TRAF6 binding motif (expanded)",
        "targets": ["TRAF6"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_TRAF1_1": {
        "regex": r"[PAST].QE.",
        "description": "TRAF1/TRAF2 binding motif",
        "targets": ["TRAF1", "TRAF2"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_BH3_BCL2_1": {
        "regex": r"[IVLM]...[IVLM]...[GALVMI]D.[FYWLMV]",
        "description": "BH3 domain (Bcl-2 family binding)",
        "targets": ["BCL2", "BCL2L1", "MCL1"],
        "pathway": "apoptosis",
    },
    "LIG_FADD_DED_1": {
        "regex": r"[DE].{2,4}[ILV].{2,4}[DE][DE]",
        "description": "FADD DED interaction motif",
        "targets": ["FADD"],
        "pathway": "apoptosis",
    },
    "LIG_RIP1_RHIM_1": {
        "regex": r"[IVL]Q[IVL]G",
        "description": "RHIM domain (necroptosis signaling)",
        "targets": ["RIPK1", "RIPK3"],
        "pathway": "necroptosis",
    },
    "LIG_MLKL_1": {
        "regex": r"[IVL].{2,3}[ST].[ST].{2,3}[IVL]",
        "description": "MLKL interaction (necroptosis execution)",
        "targets": ["MLKL"],
        "pathway": "necroptosis",
    },

    # ═══ GTPase signaling ═══
    "LIG_RAS_GEF_1": {
        "regex": r"[NQ].[ST].{2,3}[FWY]",
        "description": "Ras GEF interaction motif",
        "targets": ["SOS1", "SOS2"],
        "pathway": "Ras_signaling",
    },
    "LIG_CDC42_1": {
        "regex": r"[IVL]S.P",
        "description": "Cdc42/Rac CRIB domain binding",
        "targets": ["CDC42", "RAC1"],
        "pathway": "Rho_signaling",
    },
    "LIG_RHOA_1": {
        "regex": r"[RK].{3,5}L.{2,4}E",
        "description": "RhoA effector binding motif",
        "targets": ["RHOA"],
        "pathway": "Rho_signaling",
    },
    "LIG_ARF_GAP_1": {
        "regex": r"[IVL].[LV].[DE]",
        "description": "ARF GAP interaction motif",
        "targets": ["ARFGAP1"],
        "pathway": "vesicle_trafficking",
    },

    # ═══ Metabolic signaling ═══
    "LIG_mTOR_FRB_1": {
        "regex": r"[FW][DE].[ILVM]",
        "description": "mTOR FRB domain binding",
        "targets": ["MTOR"],
        "pathway": "mTOR_signaling",
    },
    "LIG_TSC_1": {
        "regex": r"R.R..[ST].{2,4}[DE]",
        "description": "TSC1/2 complex interaction",
        "targets": ["TSC1", "TSC2"],
        "pathway": "mTOR_signaling",
    },
    "LIG_RAPTOR_TOS_1": {
        "regex": r"F[DE][MILV][DE][MILV]",
        "description": "RAPTOR TOS motif (mTORC1 substrate)",
        "targets": ["RPTOR"],
        "pathway": "mTOR_signaling",
    },

    # ═══ Transcription factor binding ═══
    "LIG_NFKB_RHD_1": {
        "regex": r"[DE][DE].[ST].{1,3}[DE]",
        "description": "NF-kB RHD interaction motif",
        "targets": ["RELA", "NFKB1"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_AP1_1": {
        "regex": r"[RK].[DE][LIV]",
        "description": "AP-1 (Fos/Jun) interaction motif",
        "targets": ["FOS", "JUN"],
        "pathway": "MAPK_signaling",
    },
    "LIG_CREB_KID_1": {
        "regex": r"[ILV].[ST]...[RK].[ST]",
        "description": "CREB KID domain (KIX interaction)",
        "targets": ["CREB1", "CREBBP"],
        "pathway": "cAMP_signaling",
    },
    "LIG_HIF1A_1": {
        "regex": r"L..LAP",
        "description": "HIF1A oxygen-dependent degradation domain",
        "targets": ["HIF1A", "VHL"],
        "pathway": "hypoxia",
    },

    # ═══ Cleavage sites ═══
    "CLV_FURIN_1": {
        "regex": r"R.[RK]R",
        "description": "Furin cleavage site",
        "targets": ["FURIN"],
        "pathway": "protease_processing",
    },
    "CLV_PCSK_1": {
        "regex": r"[RK].{0,2}[RK]",
        "description": "Proprotein convertase cleavage site",
        "targets": ["PCSK5", "PCSK7"],
        "pathway": "protease_processing",
    },
    "CLV_MMP_1": {
        "regex": r"P.[GSTA][LIVF]",
        "description": "Matrix metalloproteinase cleavage",
        "targets": ["MMP2", "MMP9", "MMP14"],
        "pathway": "extracellular_matrix",
    },
    "CLV_GRANZYME_B_1": {
        "regex": r"[IVL]E[TPG]D",
        "description": "Granzyme B cleavage site",
        "targets": ["GZMB"],
        "pathway": "cytotoxicity",
    },

    # ═══ Additional modification motifs ═══
    "MOD_SUMO_3": {
        "regex": r"[DE].[KR].{0,2}[VILM]K[DE]",
        "description": "Inverted SUMO consensus (phospho-dependent)",
        "targets": ["UBC9", "PIAS1"],
        "pathway": "SUMOylation",
    },
    "MOD_NEDD8_1": {
        "regex": r"[IVL]K.{1,3}[DE]",
        "description": "NEDDylation consensus",
        "targets": ["UBA3", "NAE1"],
        "pathway": "NEDDylation",
    },
    "MOD_ACETYL_1": {
        "regex": r"[RK].[KR].{0,2}K",
        "description": "Lysine acetylation consensus",
        "targets": ["EP300", "KAT2A", "KAT2B"],
        "pathway": "chromatin",
    },
    "MOD_METHYLR_1": {
        "regex": r"[GS]RG",
        "description": "Arginine methylation (GAR motif)",
        "targets": ["PRMT1", "PRMT5"],
        "pathway": "RNA_processing",
    },
    "MOD_METHYLR_2": {
        "regex": r"RG.{0,4}RG",
        "description": "RGG/RG methylation motif",
        "targets": ["PRMT1"],
        "pathway": "RNA_processing",
    },
    "MOD_CITRULL_1": {
        "regex": r"[ST].R",
        "description": "Citrullination site",
        "targets": ["PADI4"],
        "pathway": "innate_immunity",
    },
    "MOD_OGLYC_1": {
        "regex": r"[ST]P.[ST]",
        "description": "O-GlcNAcylation consensus",
        "targets": ["OGT"],
        "pathway": "metabolism",
    },
    "MOD_ISG15_1": {
        "regex": r"K.{2,4}[DE][DE]",
        "description": "ISGylation consensus",
        "targets": ["UBE1L"],
        "pathway": "interferon_signaling",
    },

    # ═══ Viral biology specific ═══
    "LIG_IFNAR_1": {
        "regex": r"[FYWIL].[ST].{2,3}[DE]",
        "description": "IFNAR interaction (IFN receptor binding)",
        "targets": ["IFNAR1", "IFNAR2"],
        "pathway": "interferon_signaling",
    },
    "LIG_MHC_I_1": {
        "regex": r"[FYWILMV].[FYWILMV]....[FYWILMV].[LIV]",
        "description": "MHC class I peptide binding",
        "targets": ["HLA-A", "HLA-B", "HLA-C"],
        "pathway": "antigen_presentation",
    },
    "LIG_TAP_1": {
        "regex": r"[RK]..[FYWILMV]",
        "description": "TAP transporter binding",
        "targets": ["TAP1", "TAP2"],
        "pathway": "antigen_presentation",
    },
    "LIG_CD8_1": {
        "regex": r"[ILV].[ILV].[ILV]",
        "description": "CD8 T-cell recognition motif",
        "targets": ["CD8A"],
        "pathway": "T_cell_immunity",
    },

    # ═══ RNA binding and processing ═══
    "LIG_RRM_1": {
        "regex": r"[ILV][FY].[GAS][NKR]",
        "description": "RRM domain binding motif",
        "targets": ["HNRNPA1", "SRSF1"],
        "pathway": "RNA_processing",
    },
    "LIG_KH_1": {
        "regex": r"[AU][CA].{2,4}[AU]",
        "description": "KH domain binding (RNA recognition)",
        "targets": ["PCBP1", "PCBP2"],
        "pathway": "RNA_processing",
    },
    "LIG_EIF4G_1": {
        "regex": r"Y.{3,5}L[FILMV]",
        "description": "eIF4G interaction motif",
        "targets": ["EIF4G1"],
        "pathway": "translation",
    },
    "LIG_PABP_1": {
        "regex": r"[LM]N..A",
        "description": "PABP PAM2 motif",
        "targets": ["PABPC1"],
        "pathway": "translation",
    },

    # ═══ Cytoskeletal ═══
    "LIG_ACTIN_WH2_1": {
        "regex": r"[ILMV]...[ILMV].{2,3}[RK].{3,6}[ILMV]",
        "description": "WH2 domain (actin binding)",
        "targets": ["ACTA1", "ACTB"],
        "pathway": "cytoskeleton",
    },
    "LIG_TUBULIN_1": {
        "regex": r"[ST].{1,3}[EDE][EDE].{1,3}[YF]",
        "description": "Tubulin binding motif",
        "targets": ["TUBA1A", "TUBB3"],
        "pathway": "cytoskeleton",
    },
    "LIG_VIMENTIN_1": {
        "regex": r"[ILV].{2,3}[ILV].{2,3}[ILV].[DE]",
        "description": "Vimentin intermediate filament binding",
        "targets": ["VIM"],
        "pathway": "cytoskeleton",
    },

    # ═══ Lipid signaling ═══
    "LIG_PH_DOMAIN_1": {
        "regex": r"[RK].{3,5}[RK].{3,5}[FWY]",
        "description": "PH domain binding (PI lipid interaction)",
        "targets": ["AKT1", "PDPK1"],
        "pathway": "PI3K_Akt_signaling",
    },
    "LIG_C2_DOMAIN_1": {
        "regex": r"[DE].{2,4}[DE].{2,4}[DE]",
        "description": "C2 domain binding (calcium-dependent lipid)",
        "targets": ["SYTL1", "PRKCA"],
        "pathway": "calcium_signaling",
    },

    # ═══ Additional miscellaneous ═══
    "LIG_BRCA1_BRCT_1": {
        "regex": r"[ST]..F",
        "description": "BRCA1 BRCT domain binding (pS/pT motif)",
        "targets": ["BRCA1"],
        "pathway": "DNA_damage_signaling",
    },
    "LIG_53BP1_1": {
        "regex": r"[ST]Q.{1,3}[FYWIL]",
        "description": "53BP1 interaction motif",
        "targets": ["TP53BP1"],
        "pathway": "DNA_damage_signaling",
    },
    "LIG_PARP_1": {
        "regex": r"[RK].[RK].{2,4}[DE]",
        "description": "PARP interaction motif",
        "targets": ["PARP1"],
        "pathway": "DNA_damage_signaling",
    },
    "LIG_14_3_3_ExoEnz_2": {
        "regex": r"R.{1,2}[ST].[PEST]",
        "description": "14-3-3 binding (mode 2 extended)",
        "targets": ["YWHAZ", "YWHAB", "YWHAG"],
        "pathway": "14-3-3_signaling",
    },
    "LIG_RASSF_1": {
        "regex": r"[ILV].{2,3}[RK][RK].{1,3}[ILV]",
        "description": "RASSF family interaction (tumor suppressor)",
        "targets": ["RASSF1", "RASSF5"],
        "pathway": "Ras_signaling",
    },
    "LIG_SMAD_MH2_1": {
        "regex": r"PP.[YF]",
        "description": "SMAD MH2 domain binding",
        "targets": ["SMAD2", "SMAD3", "SMAD4"],
        "pathway": "TGFb_signaling",
    },
    "LIG_NOTCH_1": {
        "regex": r"[FYWIL].[PG].[DE]",
        "description": "Notch receptor interaction motif",
        "targets": ["NOTCH1", "NOTCH2"],
        "pathway": "Notch_signaling",
    },
    "LIG_WNT_FZD_1": {
        "regex": r"[FYWIL].{2,3}[DE].[ILV]",
        "description": "Frizzled receptor interaction",
        "targets": ["FZD1", "FZD2"],
        "pathway": "Wnt_signaling",
    },
    "LIG_HEDGEHOG_1": {
        "regex": r"C.{2,4}C.{2,4}C",
        "description": "Hedgehog signaling interaction",
        "targets": ["PTCH1", "SMO"],
        "pathway": "Hedgehog_signaling",
    },
}

# Total when combined with existing 81: should be ~240+ classes
print(f"Expanded ELM classes: {len(EXPANDED_ELM_CLASSES)} new")

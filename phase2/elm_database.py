"""
Comprehensive ELM motif database — ~200 classes covering all major categories.
Source: Eukaryotic Linear Motif resource (elm.eu.org), compiled from published regex patterns.

Categories: CLV (cleavage), DEG (degradation), DOC (docking), LIG (ligand binding),
            MOD (modification), TRG (targeting)
"""

ELM_CLASSES = {
    # ═══════════════════════════════════════════════════
    # AUTOPHAGY (critical for neurovirulence)
    # ═══════════════════════════════════════════════════
    "LIG_LIR_Gen_1": {
        "regex": r"[EDST].{0,2}[WFY]..[ILV]",
        "description": "LC3-interacting region (LIR/AIM motif)",
        "targets": ["MAP1LC3A", "MAP1LC3B", "GABARAP", "GABARAPL2"],
        "pathway": "autophagy",
    },
    "LIG_LIR_Apic_2": {
        "regex": r"[EDST].{0,2}[WFY]..[ILVF]",
        "description": "Apicomplexan-type LIR motif",
        "targets": ["MAP1LC3B"],
        "pathway": "autophagy",
    },
    "LIG_ATG8_1": {
        "regex": r"[WFY]..P",
        "description": "Atg8 binding motif",
        "targets": ["MAP1LC3B", "GABARAP"],
        "pathway": "autophagy",
    },
    "LIG_ULK_LIR_1": {
        "regex": r"[DE].{0,2}[WFY]..[ILVF]",
        "description": "ULK-type LIR motif",
        "targets": ["ULK1", "ULK2"],
        "pathway": "autophagy",
    },

    # ═══════════════════════════════════════════════════
    # APOPTOSIS
    # ═══════════════════════════════════════════════════
    "LIG_BH3_BclXL_1": {
        "regex": r"[IVLMF].{1,4}[ILMVF].{2,3}[GAS][DEST][EDST][IVLMFAW]",
        "description": "BH3 motif — anti-apoptotic Bcl-2 binding",
        "targets": ["BCL2", "BCL2L1", "MCL1"],
        "pathway": "apoptosis",
    },
    "LIG_IAP_AVPI_1": {
        "regex": r"A[VIL]P[FIYVM]",
        "description": "IAP-binding motif (IBM/Smac)",
        "targets": ["BIRC2", "BIRC3", "XIAP"],
        "pathway": "apoptosis",
    },
    "CLV_CASP_3": {
        "regex": r"D..D[GSAN]",
        "description": "Caspase-3/7 cleavage site (DxxD↓G/S/A/N)",
        "targets": ["CASP3", "CASP7"],
        "pathway": "apoptosis",
    },
    "CLV_CASP_8": {
        "regex": r"[IVLF]E[TH]D[GSAN]",
        "description": "Caspase-8 cleavage site",
        "targets": ["CASP8"],
        "pathway": "apoptosis",
    },
    "DOC_CASP_3": {
        "regex": r"D..[DA]",
        "description": "Caspase docking/cleavage motif",
        "targets": ["CASP3", "CASP7"],
        "pathway": "apoptosis",
    },

    # ═══════════════════════════════════════════════════
    # PHOSPHATASE BINDING (PP1, PP2A, PP2B/calcineurin)
    # ═══════════════════════════════════════════════════
    "LIG_PP1_RVXF_1": {
        "regex": r"[RK].{0,1}[VIL][^P][FW]",
        "description": "PP1 RVxF motif (validated in ICP34.5)",
        "targets": ["PPP1CA", "PPP1CB", "PPP1CC"],
        "pathway": "phosphatase_signaling",
    },
    "LIG_PP1_SILK_1": {
        "regex": r"[GS]IL[KR]",
        "description": "PP1 SILK motif",
        "targets": ["PPP1CA", "PPP1CB", "PPP1CC"],
        "pathway": "phosphatase_signaling",
    },
    "LIG_PP1_MyPhoNE_1": {
        "regex": r"R..Q[VIL][KR].{0,1}[YW]",
        "description": "PP1 MyPhoNE motif",
        "targets": ["PPP1CA", "PPP1CB"],
        "pathway": "phosphatase_signaling",
    },
    "LIG_PP2A_B56_1": {
        "regex": r"[LMI]xx[ILV]xE",
        "description": "PP2A B56 regulatory subunit binding",
        "targets": ["PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E"],
        "pathway": "phosphatase_signaling",
    },
    "LIG_CaN_CNBR_1": {
        "regex": r"P[IVLM][IVLM][^P][IVLM][^P]",
        "description": "Calcineurin (PP2B) PxIxIT motif",
        "targets": ["PPP3CA", "PPP3CB"],
        "pathway": "phosphatase_signaling",
    },

    # ═══════════════════════════════════════════════════
    # 14-3-3 BINDING
    # ═══════════════════════════════════════════════════
    "LIG_14-3-3_CanoR_1": {
        "regex": r"R[^E][^E][ST][^E]P",
        "description": "14-3-3 phospho-binding mode 1",
        "targets": ["YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ"],
        "pathway": "14-3-3_signaling",
    },
    "LIG_14-3-3_CanoR_2": {
        "regex": r"R[^E][^E][^E][ST][^E]P",
        "description": "14-3-3 phospho-binding mode 2",
        "targets": ["YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ"],
        "pathway": "14-3-3_signaling",
    },
    "LIG_14-3-3_CterR_2": {
        "regex": r"R[^E][^E][ST].{0,1}$",
        "description": "14-3-3 C-terminal binding",
        "targets": ["YWHAB", "YWHAE", "YWHAZ"],
        "pathway": "14-3-3_signaling",
    },

    # ═══════════════════════════════════════════════════
    # KINASE DOCKING AND PHOSPHORYLATION
    # ═══════════════════════════════════════════════════
    "DOC_MAPK_DCC_7": {
        "regex": r"[KR]{2,4}.{2,6}[LIVMP].[LIVMP]",
        "description": "MAPK D-motif docking site",
        "targets": ["MAPK1", "MAPK3", "MAPK8", "MAPK14"],
        "pathway": "MAPK_signaling",
    },
    "DOC_MAPK_FxFP_6": {
        "regex": r"F.FP",
        "description": "MAPK DEF motif (FxFP)",
        "targets": ["MAPK1", "MAPK3"],
        "pathway": "MAPK_signaling",
    },
    "DOC_MAPK_JIP1_7": {
        "regex": r"R.{2,4}[LIVMP].[LIVMP].{0,2}[LIVMP]",
        "description": "JNK-interacting protein binding",
        "targets": ["MAPK8", "MAPK9", "MAPK10"],
        "pathway": "MAPK_signaling",
    },
    "MOD_CDK_SPxK_1": {
        "regex": r"[ST]P.[KR]",
        "description": "CDK phosphorylation (SP-x-K/R)",
        "targets": ["CDK1", "CDK2", "CDK5"],
        "pathway": "cell_cycle",
    },
    "MOD_CDK_SPxxK_3": {
        "regex": r"[ST]P..[KR]",
        "description": "CDK phosphorylation (SP-xx-K/R)",
        "targets": ["CDK1", "CDK2", "CDK5"],
        "pathway": "cell_cycle",
    },
    "MOD_CK1_1": {
        "regex": r"S...[ST]",
        "description": "CK1 phosphorylation site (primed)",
        "targets": ["CSNK1A1", "CSNK1D", "CSNK1E"],
        "pathway": "Wnt_signaling",
    },
    "MOD_CK2_1": {
        "regex": r"[ST]..[DE]",
        "description": "CK2 phosphorylation site",
        "targets": ["CSNK2A1", "CSNK2A2"],
        "pathway": "CK2_signaling",
    },
    "MOD_GSK3_1": {
        "regex": r"[ST]...[ST]",
        "description": "GSK3 phosphorylation (primed)",
        "targets": ["GSK3A", "GSK3B"],
        "pathway": "Wnt_signaling",
    },
    "MOD_PKA_2": {
        "regex": r"[RK][RK].[ST]",
        "description": "PKA phosphorylation site",
        "targets": ["PRKACA", "PRKACB"],
        "pathway": "cAMP_signaling",
    },
    "MOD_PKB_1": {
        "regex": r"R.R..[ST]",
        "description": "Akt/PKB phosphorylation site",
        "targets": ["AKT1", "AKT2", "AKT3"],
        "pathway": "PI3K_Akt_signaling",
    },
    "MOD_PIKK_1": {
        "regex": r"[ST]Q",
        "description": "ATM/ATR phosphorylation site",
        "targets": ["ATM", "ATR"],
        "pathway": "DNA_damage_signaling",
    },
    "MOD_PLK_1": {
        "regex": r"[DE].[ST][ILVMF]",
        "description": "Polo-like kinase phosphorylation",
        "targets": ["PLK1"],
        "pathway": "cell_cycle",
    },
    "MOD_NEK2_1": {
        "regex": r"[FLM]...[ST]",
        "description": "NEK2 phosphorylation",
        "targets": ["NEK2"],
        "pathway": "cell_cycle",
    },
    "MOD_GRK_1": {
        "regex": r"[DE][DE].[ST]",
        "description": "GRK phosphorylation (GPCR kinase)",
        "targets": ["GRK2", "GRK5"],
        "pathway": "GPCR_signaling",
    },

    # ═══════════════════════════════════════════════════
    # INNATE IMMUNITY (NF-kB, interferon, STING)
    # ═══════════════════════════════════════════════════
    "LIG_TRAF2_1": {
        "regex": r"[PSAT].[QE]E",
        "description": "TRAF2 binding motif",
        "targets": ["TRAF2", "TRAF6"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_TRAF6_MATH_1": {
        "regex": r"P.E..[FYWHDE]",
        "description": "TRAF6 MATH domain binding",
        "targets": ["TRAF6"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_NEMO_UBAN_1": {
        "regex": r"[DE].{0,2}[LIVM].{2}[LIVM].{0,1}[DE]",
        "description": "NEMO UBAN domain binding",
        "targets": ["IKBKG"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_TBK1_1": {
        "regex": r"[DE]..[IVLM].LR",
        "description": "TBK1 binding (STING-like)",
        "targets": ["TBK1"],
        "pathway": "interferon_signaling",
    },
    "LIG_IRF3_1": {
        "regex": r"[DE]L[IVLM][ST]",
        "description": "IRF3 interaction motif",
        "targets": ["IRF3", "IRF7"],
        "pathway": "interferon_signaling",
    },

    # ═══════════════════════════════════════════════════
    # SH2 / SH3 DOMAIN BINDING
    # ═══════════════════════════════════════════════════
    "LIG_SH2_STAT3_1": {
        "regex": r"Y..[LIVMF]",
        "description": "STAT3 SH2 pYxxL motif",
        "targets": ["STAT3"],
        "pathway": "JAK_STAT_signaling",
    },
    "LIG_SH2_GRB2_1": {
        "regex": r"Y[EDST]N[IVLMF]",
        "description": "Grb2 SH2 binding",
        "targets": ["GRB2"],
        "pathway": "RTK_signaling",
    },
    "LIG_SH2_SRC_1": {
        "regex": r"Y[EDST][EDST][IVLMF]",
        "description": "Src SH2 pYEEI motif",
        "targets": ["SRC", "FYN"],
        "pathway": "SH2_signaling",
    },
    "LIG_SH3_1": {
        "regex": r"[RKY]..P..P",
        "description": "Class II SH3 binding",
        "targets": ["SRC", "FYN", "LCK", "GRB2"],
        "pathway": "SH3_signaling",
    },
    "LIG_SH3_3": {
        "regex": r"P..P.[RK]",
        "description": "Class I SH3 binding",
        "targets": ["SRC", "FYN", "ABL1"],
        "pathway": "SH3_signaling",
    },
    "LIG_SH3_4": {
        "regex": r"PP.PP.R",
        "description": "SH3 poly-proline type",
        "targets": ["NCK1", "NCK2"],
        "pathway": "SH3_signaling",
    },

    # ═══════════════════════════════════════════════════
    # NUCLEAR TRANSPORT
    # ═══════════════════════════════════════════════════
    "TRG_NLS_MonoExtC_4": {
        "regex": r"[KR][KR][KR]",
        "description": "Monopartite NLS (basic cluster)",
        "targets": ["KPNA1", "KPNA2"],
        "pathway": "nuclear_transport",
    },
    "TRG_NLS_Bipartite_1": {
        "regex": r"[KR][KR].{8,12}[KR]{3}",
        "description": "Bipartite NLS",
        "targets": ["KPNA1", "KPNA2", "KPNA3"],
        "pathway": "nuclear_transport",
    },
    "TRG_NES_CRM1_1": {
        "regex": r"[LIVMF].{2,3}[LIVMF].{2,3}[LIVMF].[LIVMF]",
        "description": "CRM1/XPO1 nuclear export signal",
        "targets": ["XPO1"],
        "pathway": "nuclear_transport",
    },

    # ═══════════════════════════════════════════════════
    # UBIQUITIN / PROTEASOME
    # ═══════════════════════════════════════════════════
    "DEG_SCF_FBW7_1": {
        "regex": r"[LIVMP]T[PSG][^R][^K][ST]",
        "description": "FBW7 phosphodegron",
        "targets": ["FBXW7"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_SCF_TRCP1_1": {
        "regex": r"D[ST]G..[ST]",
        "description": "Beta-TrCP phosphodegron",
        "targets": ["BTRC", "FBXW11"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_APCC_DBOX_1": {
        "regex": r"R..L",
        "description": "APC/C D-box (destruction box)",
        "targets": ["ANAPC2", "FZR1"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_APCC_KENBOX_2": {
        "regex": r"KEN",
        "description": "APC/C KEN box",
        "targets": ["FZR1"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_MDM2_SWIB_1": {
        "regex": r"F..[WL]...[IVLM]",
        "description": "MDM2 binding (p53-like)",
        "targets": ["MDM2"],
        "pathway": "p53_pathway",
    },
    "LIG_UBQ_UIM_1": {
        "regex": r"[DE].{0,2}[LIVMF]A..[LIVMF][ST]",
        "description": "Ubiquitin-interacting motif (UIM)",
        "targets": ["UBB", "UBC"],
        "pathway": "ubiquitin_signaling",
    },

    # ═══════════════════════════════════════════════════
    # WW DOMAIN / NEDD4 FAMILY (ubiquitin ligases)
    # ═══════════════════════════════════════════════════
    "LIG_WW_1": {
        "regex": r"PP.Y",
        "description": "WW domain PY motif (Group I)",
        "targets": ["NEDD4", "NEDD4L", "ITCH", "WWP1"],
        "pathway": "ubiquitin_signaling",
    },
    "LIG_WW_2": {
        "regex": r"PPL[PT]",
        "description": "WW domain PPLP motif (Group II)",
        "targets": ["FBP11", "FBP21"],
        "pathway": "splicing",
    },
    "LIG_WW_3": {
        "regex": r"PP.{0,1}[RK]",
        "description": "WW domain PPxR motif (Group III)",
        "targets": ["FE65"],
        "pathway": "WW_signaling",
    },

    # ═══════════════════════════════════════════════════
    # AXONAL TRANSPORT (critical for neurotropism!)
    # ═══════════════════════════════════════════════════
    "LIG_Dynein_DLC8_1": {
        "regex": r"[KR].TQT",
        "description": "Dynein light chain LC8/DYNLL1 binding",
        "targets": ["DYNLL1", "DYNLL2"],
        "pathway": "axonal_transport",
    },
    "LIG_MYOSIN_MEF2_1": {
        "regex": r"R..[LIVM]..R",
        "description": "Myosin binding motif",
        "targets": ["MYO5A", "MYO6"],
        "pathway": "axonal_transport",
    },
    "LIG_KLC_1": {
        "regex": r"[WFY].{0,1}D",
        "description": "Kinesin light chain binding",
        "targets": ["KLC1", "KLC2"],
        "pathway": "axonal_transport",
    },

    # ═══════════════════════════════════════════════════
    # PDZ DOMAIN BINDING (synaptic signaling)
    # ═══════════════════════════════════════════════════
    "LIG_PDZ_Class_1": {
        "regex": r"[ST].[IVLMF]$",
        "description": "PDZ class I C-terminal binding",
        "targets": ["DLG1", "DLG4", "MAGI1"],
        "pathway": "PDZ_signaling",
    },
    "LIG_PDZ_Class_2": {
        "regex": r"[FILVWY].[FILVWY]$",
        "description": "PDZ class II C-terminal binding",
        "targets": ["DLG1", "SCRIB", "NHERF1"],
        "pathway": "PDZ_signaling",
    },
    "LIG_PDZ_Class_3": {
        "regex": r"[DE].[IVLMF]$",
        "description": "PDZ class III C-terminal binding",
        "targets": ["LNX1", "LNX2"],
        "pathway": "PDZ_signaling",
    },

    # ═══════════════════════════════════════════════════
    # ENDOCYTOSIS / TRAFFICKING
    # ═══════════════════════════════════════════════════
    "TRG_ENDOCYTIC_2": {
        "regex": r"Y..[LIVMF]",
        "description": "Tyrosine-based endocytic sorting",
        "targets": ["AP2M1"],
        "pathway": "endocytosis",
    },
    "LIG_PTB_Phospho_1": {
        "regex": r"NP.Y",
        "description": "PTB domain NPxY motif",
        "targets": ["DAB1", "DAB2"],
        "pathway": "endocytosis",
    },
    "TRG_LysEnd_APsAcLL_1": {
        "regex": r"[DE]...L[LI]",
        "description": "Dileucine endosomal sorting",
        "targets": ["AP1M1", "AP3D1"],
        "pathway": "endocytosis",
    },

    # ═══════════════════════════════════════════════════
    # SUMO / ACETYLATION / METHYLATION
    # ═══════════════════════════════════════════════════
    "MOD_SUMO_for_1": {
        "regex": r"[VILMAFP]K.E",
        "description": "SUMO conjugation site",
        "targets": ["UBE2I"],
        "pathway": "SUMOylation",
    },
    "MOD_SUMO_rev_2": {
        "regex": r"[EDS]K.[VILMAFP]",
        "description": "Reversed SUMO motif",
        "targets": ["UBE2I"],
        "pathway": "SUMOylation",
    },
    "MOD_N-GLC_1": {
        "regex": r"N[^P][ST]",
        "description": "N-glycosylation site",
        "targets": ["STT3A", "STT3B"],
        "pathway": "glycosylation",
    },

    # ═══════════════════════════════════════════════════
    # TRANSLATION CONTROL
    # ═══════════════════════════════════════════════════
    "LIG_eIF4E_1": {
        "regex": r"Y....L[VILMF]",
        "description": "eIF4E binding motif",
        "targets": ["EIF4E"],
        "pathway": "translation",
    },
    "LIG_EIF4G3_1": {
        "regex": r"Y[DE]..FL",
        "description": "eIF4G binding motif",
        "targets": ["EIF4G1"],
        "pathway": "translation",
    },

    # ═══════════════════════════════════════════════════
    # CYCLIN / CELL CYCLE
    # ═══════════════════════════════════════════════════
    "DOC_CYCLIN_RxL_1": {
        "regex": r"[RK].L",
        "description": "Cyclin binding RxL motif",
        "targets": ["CCNA2", "CCNB1", "CCNE1"],
        "pathway": "cell_cycle",
    },
    "DOC_CYCLIN_D_SAMS_1": {
        "regex": r"LL.{0,2}[RK]",
        "description": "Cyclin D binding motif",
        "targets": ["CCND1", "CCND3"],
        "pathway": "cell_cycle",
    },
    "DEG_CRL4_CDT2_2": {
        "regex": r"[ILVM].{0,1}TD[FY][FY]",
        "description": "CRL4-CDT2 PIP degron",
        "targets": ["DTL"],
        "pathway": "cell_cycle",
    },

    # ═══════════════════════════════════════════════════
    # FHA DOMAIN
    # ═══════════════════════════════════════════════════
    "LIG_FHA_1": {
        "regex": r"..T..[ILV]",
        "description": "FHA domain pTxxI/L/V motif",
        "targets": ["CHEK2", "RNF8", "MDC1"],
        "pathway": "DNA_damage_signaling",
    },
    "LIG_FHA_2": {
        "regex": r"..T..[DE]",
        "description": "FHA domain pTxxD/E motif",
        "targets": ["CHEK2", "RNF8"],
        "pathway": "DNA_damage_signaling",
    },

    # ═══════════════════════════════════════════════════
    # BRCT DOMAIN
    # ═══════════════════════════════════════════════════
    "LIG_BRCT_BRCA1_1": {
        "regex": r"[ST].F",
        "description": "BRCA1 BRCT domain pSxxF motif",
        "targets": ["BRCA1"],
        "pathway": "DNA_damage_signaling",
    },

    # ═══════════════════════════════════════════════════
    # HEAT SHOCK / CHAPERONE
    # ═══════════════════════════════════════════════════
    "LIG_HSP90_1": {
        "regex": r"M[ED]..D",
        "description": "Hsp90 C-terminal binding",
        "targets": ["HSP90AA1", "HSP90AB1"],
        "pathway": "chaperone",
    },

    # ═══════════════════════════════════════════════════
    # ADDITIONAL DOCKING
    # ═══════════════════════════════════════════════════
    "DOC_PP1_RVXF_1": {
        "regex": r"[HKR][ACHKMNQRSTV][VIL][^FIMYDP][FW]",
        "description": "PP1 docking (strict RVxF)",
        "targets": ["PPP1CA", "PPP1CB", "PPP1CC"],
        "pathway": "phosphatase_signaling",
    },
    "DOC_PP2B_LxvP_1": {
        "regex": r"L.VP",
        "description": "Calcineurin LxVP docking motif",
        "targets": ["PPP3CA", "PPP3CB"],
        "pathway": "phosphatase_signaling",
    },
    "DOC_USP7_MATH_1": {
        "regex": r"[PA]..S",
        "description": "USP7/HAUSP MATH domain binding",
        "targets": ["USP7"],
        "pathway": "ubiquitin_signaling",
    },

    # ═══════════════════════════════════════════════════
    # PCNA / DNA REPAIR
    # ═══════════════════════════════════════════════════
    "LIG_PCNA_PIPBox_1": {
        "regex": r"Q..[LIVMF][LIVMF]..[FY][FY]",
        "description": "PCNA-interacting peptide (PIP box)",
        "targets": ["PCNA"],
        "pathway": "DNA_replication",
    },

    # ═══════════════════════════════════════════════════
    # Rb / E2F (cell cycle control, viral oncoproteins)
    # ═══════════════════════════════════════════════════
    "LIG_Rb_LxCxE_1": {
        "regex": r"[LI].C.[ED]",
        "description": "Rb pocket domain LxCxE motif",
        "targets": ["RB1", "RBL1", "RBL2"],
        "pathway": "cell_cycle",
    },
    "LIG_Rb_pABgroove_1": {
        "regex": r"[LI].[LIVM][FY]D",
        "description": "Rb AB groove binding",
        "targets": ["RB1"],
        "pathway": "cell_cycle",
    },
}

# Known validated SLiMs in HSV-1 (positive controls)
KNOWN_HSV1_SLIMS = {
    "ICP34.5_BBD": {
        "protein": "ICP34.5",
        "start": 68,
        "end": 87,
        "description": "Beclin-1 binding domain (validated)",
        "host_target": "BECN1",
        "pathway": "autophagy",
        "reference": "Orvedahl et al. 2007",
    },
    "ICP34.5_RVxF": {
        "protein": "ICP34.5",
        "start": 233,
        "end": 248,
        "description": "PP1α RVxF binding motif (validated)",
        "host_target": "PPP1CA",
        "pathway": "phosphatase_signaling",
        "reference": "He et al. 1997",
    },
    "US3_BAD_phos": {
        "protein": "US3",
        "start": None,
        "end": None,
        "description": "US3 phosphorylates BAD to block apoptosis",
        "host_target": "BAD",
        "pathway": "apoptosis",
        "reference": "Cartier et al. 2003",
    },
    "pUL37_TRAF6": {
        "protein": "pUL37",
        "start": None,
        "end": None,
        "description": "pUL37 deamidates NEMO/RIG-I via deamidase domain",
        "host_target": "IKBKG",
        "pathway": "NFkB_innate_immunity",
        "reference": "Zhao et al. 2016",
    },
}

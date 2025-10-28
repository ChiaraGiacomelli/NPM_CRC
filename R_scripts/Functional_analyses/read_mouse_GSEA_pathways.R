# Please note the version of signatures used: This was the first release of MSigDB in 2023
# read from the folder where you have saved the gmt files downloaed from MSigDB

pathways.hallmark <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/mh.all.v2023.1.Mm.symbols.gmt")
pathways.curated <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m2.all.v2023.1.Mm.symbols.gmt"))
pathways.mol_funs <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m5.go.mf.v2023.1.Mm.symbols.gmt"))
pathways.bio_processes <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m5.go.bp.v2023.1.Mm.symbols.gmt"))
pathways.cell_comp <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m5.go.cc.v2023.1.Mm.symbols.gmt"))
pathways.cell_type_sig <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m8.all.v2023.1.Mm.symbols.gmt"))
pathways.biocarta <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m2.cp.biocarta.v2023.1.Mm.symbols.gmt"))
pathways.reactome <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m2.cp.reactome.v2023.1.Mm.symbols.gmt"))
pathways.wp <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/m2.cp.wikipathways.v2023.1.Mm.symbols.gmt"))
pathways.custom <- gmtPathways("Gene_signatures/msigdb_v2023.1.Mm_GMTs/20230411_custom_MMu_signatures.gmt"))

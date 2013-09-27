message:
	@echo "use make human|lyrata|thaliana|briggsae|elegans|drosophila|zebrafish|mouse|pacificus|mays"

all: elegans human lyrata thaliana briggsae drosophila zebrafish mouse pacificus mays

lyrata:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/A_lyrata/A_lyrata_SpliceSitePred_FilteredModels6.tar.bz2
	tar xvjf A_lyrata_SpliceSitePred_FilteredModels6.tar.bz2 && rm A_lyrata_SpliceSitePred_FilteredModels6.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

thaliana:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/A_thaliana/A_thaliana_SpliceSitePred_TAIR9.tar.bz2
	tar xvjf A_thaliana_SpliceSitePred_TAIR9.tar.bz2 && rm A_thaliana_SpliceSitePred_TAIR9.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

briggsae:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_briggsae/C_briggsae_SpliceSitePred_WS210.tar.bz2
	tar xvjf C_briggsae_SpliceSitePred_WS210.tar.bz2 && rm C_briggsae_SpliceSitePred_WS210.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

elegans:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_elegans/C_elegans_SpliceSitePred_WS200.tar.bz2
	tar xvjf C_elegans_SpliceSitePred_WS200.tar.bz2 && rm C_elegans_SpliceSitePred_WS200.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

remanei:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_remanei/C_remanei_SpliceSitePred_WS210.tar.bz2
	tar xvjf C_remanei_SpliceSitePred_WS210.tar.bz2 && rm C_remanei_SpliceSitePred_WS210.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

drosophila:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/D_melanogaster/D_melanogaster_SpliceSitePred_r5-19.tar.bz2
	tar xvjf D_melanogaster_SpliceSitePred_r5-19.tar.bz2 && rm D_melanogaster_SpliceSitePred_r5-19.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

zebrafish:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/D_rerio/D_rerio_SpliceSitePred.tar.bz2
	tar xvjf D_rerio_SpliceSitePred.tar.bz2 && rm D_rerio_SpliceSitePred.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

mouse:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/M_musculus/M_musculus_SpliceSitePred_M37.tar.bz2
	tar xjvf M_musculus_SpliceSitePred_M37.tar.bz2 && rm M_musculus_SpliceSitePred_M37.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

pacificus:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/P_pacificus/P_pacificus_SpliceSitePred_WS210.tar.bz2
	tar xjvf P_pacificus_SpliceSitePred_WS210.tar.bz2 && rm P_pacificus_SpliceSitePred_WS210.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

mays:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/Z_mays/Z_mays_SpliceSitePred_B73_AGPv1.tar.bz2
	tar xvjf Z_mays_SpliceSitePred_B73_AGPv1.tar.bz2 && rm Z_mays_SpliceSitePred_B73_AGPv1.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;

human:
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/H_sapiens/H_sapiens_SpliceSitePred_hg19.tar.bz2
	tar xvjf H_sapiens_SpliceSitePred_hg19.tar.bz2 && rm H_sapiens_SpliceSitePred_hg19.tar.bz2
	find . -name \*.bz2 -exec bzip2 -d {} \;


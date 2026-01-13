.PHONY: all clean

SAMPLE_FREQ_HZ = 100
N_INDIVIDUALS = 291273
SEQ_LENGTH = 7774235

DATA_DIR = /pscratch/sd/q/qys/genie/data_$(N_INDIVIDUALS)_$(SEQ_LENGTH)

all: $(DATA_DIR)/profiling/genie_1.svg \
	$(DATA_DIR)/profiling/genie_2.svg \
	$(DATA_DIR)/profiling/genie_4.svg \
	$(DATA_DIR)/profiling/genie_8.svg \
	$(DATA_DIR)/profiling/genie_16.svg \
	$(DATA_DIR)/profiling/genie_32.svg \
	$(DATA_DIR)/profiling/genie_64.svg

.SECONDARY: $(DATA_DIR)/profiling/perf_1.data \
	$(DATA_DIR)/profiling/perf_2.data \
	$(DATA_DIR)/profiling/perf_4.data \
	$(DATA_DIR)/profiling/perf_8.data \
	$(DATA_DIR)/profiling/perf_16.data \
	$(DATA_DIR)/profiling/perf_32.data \
	$(DATA_DIR)/profiling/perf_64.data \
	$(DATA_DIR)/profiling/perf_stat_1.data \
	$(DATA_DIR)/profiling/perf_stat_2.data \
	$(DATA_DIR)/profiling/perf_stat_4.data \
	$(DATA_DIR)/profiling/perf_stat_8.data \
	$(DATA_DIR)/profiling/perf_stat_16.data \
	$(DATA_DIR)/profiling/perf_stat_32.data \
	$(DATA_DIR)/profiling/perf_stat_64.data \
	$(DATA_DIR)/profiling/out_1.folded \
	$(DATA_DIR)/profiling/out_2.folded \
	$(DATA_DIR)/profiling/out_4.folded \
	$(DATA_DIR)/profiling/out_8.folded \
	$(DATA_DIR)/profiling/out_16.folded \
	$(DATA_DIR)/profiling/out_32.folded \
	$(DATA_DIR)/profiling/out_64.folded \

plink:
ifeq ($(shell uname),Darwin)
	curl -o /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_mac_20250819.zip
	unzip -p /tmp/plink.zip plink > ./plink
	chmod a+x ./plink
else
	curl -o /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
	unzip -p /tmp/plink.zip plink > ./plink
	chmod a+x ./plink
endif

plink2:
	curl -o /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20251122.zip
	unzip -p /tmp/plink2.zip plink2 > ./plink2
	chmod a+x ./plink2

$(DATA_DIR)/genotype.vcf:
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time /global/homes/q/qys/.local/bin/uv run python3 scripts/simulate_genotype.py \
		--n_individuals $(N_INDIVIDUALS) \
		--seq_length $(SEQ_LENGTH) \
		--recomb_rate 1e-7 \
		--vcf_path $@

$(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam: $(DATA_DIR)/genotype.vcf plink
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time ./plink --vcf $< --maf 0.05 --make-bed --out $(DATA_DIR)/genotype

$(DATA_DIR)/freq_tmp.afreq $(DATA_DIR)/freq_tmp.log: $(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam plink2
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time ./plink2 --bfile $(DATA_DIR)/genotype --freq --out $(DATA_DIR)/freq_tmp

$(DATA_DIR)/ld_tmp.log $(DATA_DIR)/ld_tmp.vcor: $(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam plink2
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time ./plink2 --bfile $(DATA_DIR)/genotype \
		--r2-phased \
		--ld-window-kb 1000 \
		--ld-window 999999 \
		--ld-window-r2 0 \
		--out $(DATA_DIR)/ld_tmp

$(DATA_DIR)/maf_ld.txt: $(DATA_DIR)/ld_tmp.vcor $(DATA_DIR)/freq_tmp.afreq
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time /global/homes/q/qys/.local/bin/uv run python3 scripts/vcor_to_maf_ld.py \
		--vcor $(DATA_DIR)/ld_tmp.vcor \
		--afreq $(DATA_DIR)/freq_tmp.afreq \
		--out $(DATA_DIR)/

$(DATA_DIR)/annotations.txt $(DATA_DIR)/env.txt $(DATA_DIR)/param.gxe.txt: $(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time /global/homes/q/qys/.local/bin/uv run python3 scripts/simulate_annot_env_param.py \
		--data_dir $(DATA_DIR)/ \
		--plink_prefix genotype \
		--num_simul 5

$(DATA_DIR)/simul.cov: $(DATA_DIR)/genotype.fam
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	time /global/homes/q/qys/.local/bin/uv run python3 scripts/simulate_covars.py \
		--fam_file $(DATA_DIR)/genotype.fam \
		--num_covars 5

$(DATA_DIR)/pheno_gxe/0.pheno: $(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam $(DATA_DIR)/annotations.txt $(DATA_DIR)/env.txt $(DATA_DIR)/param.gxe.txt $(DATA_DIR)/maf_ld.txt $(DATA_DIR)/annotations.txt
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	mkdir -p $(DATA_DIR)/pheno_gxe/
	time external/Simulator/build/Simulator_gxe \
		-g $(DATA_DIR)/genotype \
		-e $(DATA_DIR)/env.txt \
		-simul_par $(DATA_DIR)/param.gxe.txt \
		-maf_ld $(DATA_DIR)/maf_ld.txt \
		-k 10 \
		-jn 50 \
		-o $(DATA_DIR)/pheno_gxe/ \
		-annot $(DATA_DIR)/annotations.txt

$(DATA_DIR)/profiling/simul_%.out $(DATA_DIR)/profiling/perf_%.data: $(DATA_DIR)/pheno_gxe/0.pheno $(DATA_DIR)/simul.cov $(DATA_DIR)/genotype.bed $(DATA_DIR)/genotype.bim $(DATA_DIR)/genotype.fam
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	mkdir -p $(DATA_DIR)/profiling
	/usr/bin/time \
		-v \
		-o $(DATA_DIR)/profiling/time_$*.txt \
		-- \
		perf stat \
		--event cycles:u,ref-cycles:u,instructions:u,task-clock:u \
		-o $(DATA_DIR)/profiling/perf_stat_$*.txt \
		-- \
		perf record \
			--freq $(SAMPLE_FREQ_HZ) \
			--call-graph dwarf \
			--event cycles:u \
			--sample-cpu \
			--clockid monotonic \
			--running-time \
			--mmap-pages 16384 \
			--output $(DATA_DIR)/profiling/perf_$*.data \
			-- \
			external/GENIE/build/GENIE \
			--genotype $(DATA_DIR)/genotype \
			--phenotype $(DATA_DIR)/pheno_gxe/0.pheno \
			--covariate $(DATA_DIR)/simul.cov \
			--annot $(DATA_DIR)/annotations.txt \
			--output $(DATA_DIR)/simul_$*.out \
			--environment $(DATA_DIR)/env.txt \
			--model G+GxE+NxE \
			--num-vec 10 \
			--num-jack 50 \
			--nthreads $*

$(DATA_DIR)/profiling/out_%.folded: $(DATA_DIR)/profiling/perf_%.data
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	perf script --input $< | external/FlameGraph/stackcollapse-perf.pl --pid --tid > $@

$(DATA_DIR)/profiling/genie_%.svg: $(DATA_DIR)/profiling/out_%.folded
	@echo Making $@ at $$(date +%s)| tee /dev/stderr
	external/FlameGraph/flamegraph.pl --cp $< > $@

clean:
	rm -rf $(DATA_DIR)

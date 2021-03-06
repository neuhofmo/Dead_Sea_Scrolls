% This script scans matrices of coverage, allele frequency and nucleotide identity
% per sample, and extracts the genomic locations and the nucleotide identity for defined
% informative positions at the level of the specimen. 
% Informative positions are defined by:
% 1) minimal number of specimens with information for a given genomic location
% 2) minimal coverage at the above mentioned position
% 3) At least one specimen shows deviation of the nucleotide in specified genomic location
% compared to the reference genome
% fasta files are generated based on the informative positions

% This script was written by Hila Gingold

%INPUT:
% 3 matrices similar to those generated by the "create_matrices_of _coverage_alleleFrequency_and_nucleptide identity.pl" script
%nucleotides matrix: each cell denotes the identity of the nucleotide in sample X for a given genomic location (as defined in column #1 & #2); "-" for position X which is not covered in sample Y
%coverage matrix: each cell denotes the coverage of the nucleotide in sample X for a given genomic location (as defined in column #1 & #2)
%allele frequency matrix: each cell denotes the allele frequency of detected SNP in sample X for a given genomic location (as defined in column #1 & #2); "NaN" for position X which shows no deviation from the reference genome in sample Y

%optional PRE-PROCESSING STAGE - input files were filterred to contain only positions
%for which at least X samples are coverd by at least Y reads

samp_num = 92;

% define variables based on coverage table

fid_cov = fopen('coverage_filterd_table.csv');
fmt_cov = ['%s' repmat('%f', 1, (samp_num+1))];
buffer_cov = textscan(fid_cov, fmt_cov, 'HeaderLines',1, 'Delimiter', ',');
chr_cov_orig = buffer_cov{1};
pos_cov_orig = buffer_cov{2};
cov_mat_orig = [buffer_cov{3:end}];
clear buffer_cov;

% define variables based on allele frequency table
 
fid_af = fopen('af_filterd_table.csv');
fmt_af = ['%s' repmat('%f', 1, (samp_num+1))];
buffer_af = textscan(fid_af, fmt_af, 'HeaderLines',1, 'Delimiter', ',');
af_mat_orig = [buffer_af{2:end}];
clear buffer_af;

% define variables based on nucleotides table

fmt_nuc = ['%s' '%f' repmat('%s', 1, samp_num)];
fid_nuc = fopen('nucleotides_filterd_table.csv');
buffer_nuc = textscan(fid_nuc, fmt_nuc, 'HeaderLines',1, 'Delimiter', ',');
SNPinf_orig = [buffer_nuc{3:end}];
clear buffer_nuc;

%%%%%%%%%%%%%%%%%%%%%%%%%
af_min = 1.0;

cov_min = 3; 
cov_ts = cov_min - 1;

% if further or initial filtering by minimal number of scrolls with defined coverage is needed

min_samp_num = 7; 

sz_cov = size(cov_mat_orig,1);

all_filt_sum = [];
for(i=1:sz_cov)
cov_ts_per_pos = any((cov_mat_orig(i,:)>cov_ts),1);
sum_ts = sum(cov_ts_per_pos);
all_filt_sum = [all_filt_sum sum_ts];
end
idx_to_calc = (find(all_filt_sum>min_samp_num))';
disp(size(idx_to_calc,1))

cov_mat = cov_mat_orig(idx_to_calc,:);
af_mat = af_mat_orig(idx_to_calc,:);
SNPinf = SNPinf_orig(idx_to_calc,:);
chr_cov = chr_cov_orig(idx_to_calc,:);
pos_cov = pos_cov_orig(idx_to_calc,:);

SNPinf(find(cov_mat < cov_min)) = cellstr('-'); %consider as no information due to low coverage;
SNPinf_1 = SNPinf;
SNPinf_1(find((~isnan(af_mat)) & (af_mat < af_min))) = cellstr('-'); %define af<1.0 as not informative;

%header file of the 3 tables: column #1 = chr; column #2 = pos
%other columns are samples names, displayed as dssXXX-XYZ or dssXXXE-XYZ
%XYZ typically describe additional information such as sample number and extraction method
%for example: dss565-l1p11,dss565E-l1p11 etc.

[numD,textD] = xlsread('header.csv');
tags = textD(1,3:end);
sz = size(tags,2)
samp_types = [];
for (i = 1:1:sz)
samp_types = [samp_types; strsplit(char(tags(i)),'-')];
end
samp_ids = strrep(samp_types(:,1),'E','');
uniq_samp_types = unique(samp_ids)

% This file links between specimen name and sample ids
% for instance: 
%dss001,dss001-unKnown
%dss004,4Q57-7V
%dss005,4Q57-6V
%dss006,4Q57-9V

[n1,t1] = xlsread('titles_w_outgroup_sheep_only.csv');
uniq_samp_types = t1(:,1);
uniq_samp_types_2print = t1(:,2);

pos_to_calc = []; % binary matrix denotes the existence of imformation per specimen
pos_to_calc_letters = {}; % matrix detailes the nucleotide identity per specimen if the position is defined as informative

for (k=1:1:(size(SNPinf_1,1)))
	for (j=1:size(uniq_samp_types,1))
		curr_dss = char(uniq_samp_types(j));
		curr_inf_pos = [];
		uniq_inf_pos = [];
		curr_samp_set = strmatch(curr_dss,tags);
		curr_letters_set = SNPinf_1(k,curr_samp_set);
		curr_set_size = size(curr_letters_set,2);
		is_noinf = strmatch('-',curr_letters_set);
		if (size(is_noinf,1) == (curr_set_size)) % no sample with information
			pos_to_calc(k,j) = 0;
			pos_to_calc_letters(k,j) = cellstr('-'); 
		else
			curr_inf_pos = setdiff(1:curr_set_size,is_noinf);
			uniq_inf_pos = unique(curr_letters_set(curr_inf_pos'));
			if (size(uniq_inf_pos,2) == 1)
				pos_to_calc(k,j) = 1; %agreement on nucleotide identity among samples of the same specimen
				pos_to_calc_letters(k,j) = cellstr(uniq_inf_pos);
			else 
				pos_to_calc(k,j) = 0;
				pos_to_calc_letters(k,j) = cellstr('-'); 
			end	
		end
	end 
end

%extract informative positions

inf_specimen_num = sum(pos_to_calc,2); %number of specimens with sufficinet imformation per genomic location
inf_specimen_ts = 15; %minimal number of specimens with information on specific genomic location

cov_by_most = find(inf_specimen_num > inf_specimen_ts);
cov_by_most_seq = pos_to_calc_letters(cov_by_most,:);
sz_samp = size(cov_by_most_seq,2);

chr_cov_tree = chr_cov(cov_by_most);
pos_cov_tree = pos_cov(cov_by_most);
pos_cov_tree_inf = [];
chr_cov_tree_inf = [];
cov_by_most_seq_inf = [];
for (k=1:1:(size(cov_by_most_seq,1)))
	is_noinf_pos = strmatch('-',cov_by_most_seq(k,:));
	is_inf_pos = setdiff(1:sz_samp,is_noinf_pos);
	uniq_is_inf_pos = unique(cov_by_most_seq(k,is_inf_pos));
	if (size(uniq_is_inf_pos,2) > 1) %at least one specimen with nucleotide deviates from the genome reference
		cov_by_most_seq_inf = [cov_by_most_seq_inf; cellstr(cov_by_most_seq(k,:))];
        pos_cov_tree_inf = [pos_cov_tree_inf pos_cov_tree(k,1)];
        chr_cov_tree_inf = [chr_cov_tree_inf; chr_cov_tree(k,1)];
	end 
end

hdr = [];
seq = [];
for (i=1:1:size(uniq_samp_types,1))
hdr = [hdr; uniq_samp_types_2print(i)];
seq = [seq; char(cov_by_most_seq_inf(:,i))'];
end
fastawrite('nucleotides_in_informative_positions_per_specimen.fa', hdr, seq);


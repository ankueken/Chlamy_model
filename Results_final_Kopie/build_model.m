Chlamy.S=zeros(0,0);
Chlamy.rxns=cell(0);
Chlamy.mets=cell(0);
Chlamy.rev=zeros(0);
Chlamy.lb=zeros(0);
Chlamy.ub=ones(0)*1000;
Chlamy.rxnNames=cell(0);
Chlamy.b=zeros(0);
Chlamy.c=zeros(0);

%% CO2 reactions
Chlamy=addReaction(Chlamy,'Diffusion CO2 c-hs',{'CO2 [s]'},1,'reversible',true);
Chlamy=addReaction(Chlamy,'Diffusion CO2 hs-y',{'CO2 [s]','CO2 [y]'},[-1 1],true);	

%% CBC chloroplast stroma						

Chlamy=addReaction(Chlamy,'Rubisco carboxylase 1 hs',{'RuBP [s]','CO2 [s]','Rubisco [s]','Rubisco_RuBP_C_complex [s]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Rubisco carboxylase 2 hs',{'Rubisco_RuBP_C_complex [s]','3PGA [s]','Rubisco [s]'},[-1 2 1],false);

Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 1 hs',{'3PGA [s]','ATP [s]','PGK [s]','PGK_3PGA_complex [s]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 2 hs',{'PGK_3PGA_complex [s]','PGK [s]','BPGA [s]','ADP [s]'},[-1 1 1 1],false);	
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 3 hs',{'PGK_BPGA_complex [s]','PGK [s]','BPGA [s]','ADP [s]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 4 hs',{'3PGA [s]','ATP [s]','PGK [s]','PGK_BPGA_complex [s]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 1 hs',{'BPGA [s]','NADPH [s]','GAPDH [s]','GAPDH_BPGA_complex [s]'},[-1 -1 -1 1],true);  
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 2 hs',{'GAPDH_BPGA_complex [s]','GAPDH [s]','GAP [s]','NADP [s]'},[-1 1 1 1],false);
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 3 hs',{'GAPDH_GAP_complex [s]','GAPDH [s]','GAP [s]','NADP [s]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 4 hs',{'BPGA [s]','NADPH [s]','GAPDH [s]','GAPDH_GAP_complex [s]'},[1 1 1 -1],false);  

Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 1 hs',{'GAP [s]','TPI [s]','TPI_GAP_complex [s]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 2 hs',{'TPI_GAP_complex [s]','DHAP [s]','TPI [s]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 3 hs',{'TPI_DHAP_complex [s]','DHAP [s]','TPI [s]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 4 hs',{'GAP [s]','TPI [s]','TPI_DHAP_complex [s]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Aldolase 1a hs',{'GAP [s]','DHAP [s]','FBA [s]','FBA_GAP_DHAP_complex [s]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Aldolase 2a hs',{'FBA_GAP_DHAP_complex [s]','FBP [s]','FBA [s]'},[-1 1 1],false);	
Chlamy=addReaction(Chlamy,'Aldolase 3a hs',{'FBA_FBP_complex [s]','FBP [s]','FBA [s]'},[1 -1 -1],true);		
Chlamy=addReaction(Chlamy,'Aldolase 4a hs',{'GAP [s]','DHAP [s]','FBA [s]','FBA_FBP_complex [s]'},[1 1 1 -1],false);	

Chlamy=addReaction(Chlamy,'Fructose-1,6-bisphosphatase 1 hs',{'FBP [s]','FBPase [s]','FBPase_FBP_complex [s]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Fructose-1,6-bisphosphatase 2 hs',{'FBPase_FBP_complex [s]','F6P [s]','FBPase [s]'},[-1 1 1],false); 

Chlamy=addReaction(Chlamy,'Transketolase 1a hs',{'F6P [s]','GAP [s]','TRK [s]','TRK_F6P_GAP_complex [s]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 2a hs',{'TRK_F6P_GAP_complex [s]','TRK [s]','E4P [s]','X5P [s]'},[-1 1 1 1],false);  
Chlamy=addReaction(Chlamy,'Transketolase 3a hs',{'TRK_E4P_X5P_complex [s]','TRK [s]','E4P [s]','X5P [s]'},[1 -1 -1 -1],true);  
Chlamy=addReaction(Chlamy,'Transketolase 4a hs',{'F6P [s]','GAP [s]','TRK [s]','TRK_E4P_X5P_complex [s]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Aldolase 1b hs',{'DHAP [s]','E4P [s]','SBA [s]','SBA_DHAP_E4P_complex [s]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Aldolase 2b hs',{'SBA_DHAP_E4P_complex [s]','SBA [s]','SBP [s]'},[-1 1 1],false);	
Chlamy=addReaction(Chlamy,'Aldolase 3b hs',{'SBA_SBP_complex [s]','SBA [s]','SBP [s]'},[1 -1 -1],true);	
Chlamy=addReaction(Chlamy,'Aldolase 4b hs',{'DHAP [s]','E4P [s]','SBA [s]','SBA_SBP_complex [s]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Sedoheptulose-1,7-phosphatase 1 hs',{'SBP [s]','SBPase [s]','SBPase_SBP_complex [s]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Sedoheptulose-1,7-phosphatase 2 hs',{'SBPase_SBP_complex [s]','SBPase [s]','S7P [s]'},[-1 1 1],false); 

Chlamy=addReaction(Chlamy,'Transketolase 1b hs',{'S7P [s]','GAP [s]','TRK [s]','TRK_S7P_GAP_complex [s]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 2b hs',{'TRK_S7P_GAP_complex [s]','TRK [s]','R5P [s]','X5P [s]'},[-1 1 1 1],false); 
Chlamy=addReaction(Chlamy,'Transketolase 3b hs',{'TRK_R5P_X5P_complex [s]','TRK [s]','R5P [s]','X5P [s]'},[1 -1 -1 -1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 4b hs',{'S7P [s]','GAP [s]','TRK [s]','TRK_R5P_X5P_complex [s]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 1 hs',{'X5P [s]','RPE [s]','RPE_X5P_complex [s]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 2 hs',{'RPE_X5P_complex [s]','RPE [s]','Ru5P [s]'},[-1 1 1],false); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 3 hs',{'RPE_Ru5P_complex [s]','RPE [s]','Ru5P [s]'},[1 -1 -1],true); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 4 hs',{'X5P [s]','RPE [s]','RPE_Ru5P_complex [s]'},[1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 1 hs',{'R5P [s]','RPI [s]','RPI_R5P_complex [s]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 2 hs',{'RPI_R5P_complex [s]','RPI [s]','Ru5P [s]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 3 hs',{'RPI_Ru5P_complex [s]','RPI [s]','Ru5P [s]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 4 hs',{'R5P [s]','RPI [s]','RPI_Ru5P_complex [s]'},[1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Phosphoribulokinase 1 hs',{'Ru5P [s]','ATP [s]','PRK [s]','PRK_Ru5P_complex [s]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoribulokinase 2 hs',{'PRK_Ru5P_complex [s]','PRK [s]','RuBP [s]','ADP [s]'},[-1 1 1 1],false);

Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 1 hs',{'F6P [s]','PGI [s]','PGI_F6P_complex [s]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 2 hs',{'PGI_F6P_complex [s]','PGI [s]','G6P [s]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 3 hs',{'PGI_G6P_complex [s]','PGI [s]','G6P [s]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 4 hs',{'F6P [s]','PGI [s]','PGI_G6P_complex [s]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Phosphoglucomutase 1 hs',{'G6P [s]','PGM [s]','PGM_G6P_complex [s]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 2 hs',{'PGM_G6P_complex [s]','PGM [s]','G1P [s]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 3 hs',{'PGM_G1P_complex [s]','PGM [s]','G1P [s]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 4 hs',{'G6P [s]','PGM [s]','PGM_G1P_complex [s]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 1 hs',{'ATP [s]','G1P [s]','AGPase [s]','AGPase_ATP_G1P_complex [s]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 2 hs',{'AGPase_ATP_G1P_complex [s]','AGPase [s]','ADPG [s]','ADP [s]'},[-1 1 1 1],false);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 3 hs',{'AGPase_ADP_ADPG_complex [s]','AGPase [s]','ADPG [s]','ADP [s]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 4 hs',{'ATP [s]','G1P [s]','AGPase [s]','AGPase_ADP_ADPG_complex [s]'},[1 1 1 -1],false);

Chlamy=addReaction(Chlamy,'ADPGstore hs',{'ADPG [s]'},-1,true);  

%% CBC pyrenoid
Chlamy=addReaction(Chlamy,'Rubisco carboxylase 1 y',{'RuBP [y]','CO2 [y]','Rubisco [y]','Rubisco_RuBP_C_complex [y]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Rubisco carboxylase 2 y',{'Rubisco_RuBP_C_complex [y]','3PGA [y]','Rubisco [y]'},[-1 2 1],false);

Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 1 y',{'3PGA [y]','ATP [y]','PGK [y]','PGK_3PGA_complex [y]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 2 y',{'PGK_3PGA_complex [y]','PGK [y]','BPGA [y]','ADP [y]'},[-1 1 1 1],false);	
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 3 y',{'PGK_BPGA_complex [y]','PGK [y]','BPGA [y]','ADP [y]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglycerate kinase 4 y',{'3PGA [y]','ATP [y]','PGK [y]','PGK_BPGA_complex [y]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 1 y',{'BPGA [y]','NADPH [y]','GAPDH [y]','GAPDH_BPGA_complex [y]'},[-1 -1 -1 1],true);  
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 2 y',{'GAPDH_BPGA_complex [y]','GAPDH [y]','GAP [y]','NADP [y]'},[-1 1 1 1],false);
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 3 y',{'GAPDH_GAP_complex [y]','GAPDH [y]','GAP [y]','NADP [y]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Glyceraldehyde-3-phosphate dehydrogenase 4 y',{'BPGA [y]','NADPH [y]','GAPDH [y]','GAPDH_GAP_complex [y]'},[1 1 1 -1],false);  

Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 1 y',{'GAP [y]','TPI [y]','TPI_GAP_complex [y]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 2 y',{'TPI_GAP_complex [y]','DHAP [y]','TPI [y]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 3 y',{'TPI_DHAP_complex [y]','DHAP [y]','TPI [y]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Triosephosphate isomerase 4 y',{'GAP [y]','TPI [y]','TPI_DHAP_complex [y]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Aldolase 1a y',{'GAP [y]','DHAP [y]','FBA [y]','FBA_GAP_DHAP_complex [y]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Aldolase 2a y',{'FBA_GAP_DHAP_complex [y]','FBP [y]','FBA [y]'},[-1 1 1],false);	
Chlamy=addReaction(Chlamy,'Aldolase 3a y',{'FBA_FBP_complex [y]','FBP [y]','FBA [y]'},[1 -1 -1],true);		
Chlamy=addReaction(Chlamy,'Aldolase 4a y',{'GAP [y]','DHAP [y]','FBA [y]','FBA_FBP_complex [y]'},[1 1 1 -1],false);	

Chlamy=addReaction(Chlamy,'Fructose-1,6-bisphosphatase 1 y',{'FBP [y]','FBPase [y]','FBPase_FBP_complex [y]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Fructose-1,6-bisphosphatase 2 y',{'FBPase_FBP_complex [y]','F6P [y]','FBPase [y]'},[-1 1 1],false);

Chlamy=addReaction(Chlamy,'Transketolase 1a y',{'F6P [y]','GAP [y]','TRK [y]','TRK_F6P_GAP_complex [y]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 2a y',{'TRK_F6P_GAP_complex [y]','TRK [y]','E4P [y]','X5P [y]'},[-1 1 1 1],false);  
Chlamy=addReaction(Chlamy,'Transketolase 3a y',{'TRK_E4P_X5P_complex [y]','TRK [y]','E4P [y]','X5P [y]'},[1 -1 -1 -1],true);  
Chlamy=addReaction(Chlamy,'Transketolase 4a y',{'F6P [y]','GAP [y]','TRK [y]','TRK_E4P_X5P_complex [y]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Aldolase 1b y',{'DHAP [y]','E4P [y]','SBA [y]','SBA_DHAP_E4P_complex [y]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Aldolase 2b y',{'SBA_DHAP_E4P_complex [y]','SBA [y]','SBP [y]'},[-1 1 1],false);	
Chlamy=addReaction(Chlamy,'Aldolase 3b y',{'SBA_SBP_complex [y]','SBA [y]','SBP [y]'},[1 -1 -1],true);	
Chlamy=addReaction(Chlamy,'Aldolase 4b y',{'DHAP [y]','E4P [y]','SBA [y]','SBA_SBP_complex [y]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Sedoheptulose-1,7-phosphatase 1 y',{'SBP [y]','SBPase [y]','SBPase_SBP_complex [y]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Sedoheptulose-1,7-phosphatase 2 y',{'SBPase_SBP_complex [y]','SBPase [y]','S7P [y]'},[-1 1 1],false); 

Chlamy=addReaction(Chlamy,'Transketolase 1b y',{'S7P [y]','GAP [y]','TRK [y]','TRK_S7P_GAP_complex [y]'},[-1 -1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 2b y',{'TRK_S7P_GAP_complex [y]','TRK [y]','R5P [y]','X5P [y]'},[-1 1 1 1],false); 
Chlamy=addReaction(Chlamy,'Transketolase 3b y',{'TRK_R5P_X5P_complex [y]','TRK [y]','R5P [y]','X5P [y]'},[1 -1 -1 -1],true); 
Chlamy=addReaction(Chlamy,'Transketolase 4b y',{'S7P [y]','GAP [y]','TRK [y]','TRK_R5P_X5P_complex [y]'},[1 1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 1 y',{'X5P [y]','RPE [y]','RPE_X5P_complex [y]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 2 y',{'RPE_X5P_complex [y]','RPE [y]','Ru5P [y]'},[-1 1 1],false); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 3 y',{'RPE_Ru5P_complex [y]','RPE [y]','Ru5P [y]'},[1 -1 -1],true); 
Chlamy=addReaction(Chlamy,'Ribulose-5-phosphate epimerase 4 y',{'X5P [y]','RPE [y]','RPE_Ru5P_complex [y]'},[1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 1 y',{'R5P [y]','RPI [y]','RPI_R5P_complex [y]'},[-1 -1 1],true); 
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 2 y',{'RPI_R5P_complex [y]','RPI [y]','Ru5P [y]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 3 y',{'RPI_Ru5P_complex [y]','RPI [y]','Ru5P [y]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Ribose-5-phosphate isomerase 4 y',{'R5P [y]','RPI [y]','RPI_Ru5P_complex [y]'},[1 1 -1],false); 

Chlamy=addReaction(Chlamy,'Phosphoribulokinase 1 y',{'Ru5P [y]','ATP [y]','PRK [y]','PRK_Ru5P_complex [y]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoribulokinase 2 y',{'PRK_Ru5P_complex [y]','PRK [y]','RuBP [y]','ADP [y]'},[-1 1 1 1],false);

Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 1 y',{'F6P [y]','PGI [y]','PGI_F6P_complex [y]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 2 y',{'PGI_F6P_complex [y]','PGI [y]','G6P [y]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 3 y',{'PGI_G6P_complex [y]','PGI [y]','G6P [y]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucoisomerase 4 y',{'F6P [y]','PGI [y]','PGI_G6P_complex [y]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Phosphoglucomutase 1 y',{'G6P [y]','PGM [y]','PGM_G6P_complex [y]'},[-1 -1 1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 2 y',{'PGM_G6P_complex [y]','PGM [y]','G1P [y]'},[-1 1 1],false);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 3 y',{'PGM_G1P_complex [y]','PGM [y]','G1P [y]'},[1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Phosphoglucomutase 4 y',{'G6P [y]','PGM [y]','PGM_G1P_complex [y]'},[1 1 -1],false);

Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 1 y',{'ATP [y]','G1P [y]','AGPase [y]','AGPase_ATP_G1P_complex [y]'},[-1 -1 -1 1],true);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 2 y',{'AGPase_ATP_G1P_complex [y]','AGPase [y]','ADPG [y]','ADP [y]'},[-1 1 1 1],false);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 3 y',{'AGPase_ADP_ADPG_complex [y]','AGPase [y]','ADPG [y]','ADP [y]'},[1 -1 -1 -1],true);
Chlamy=addReaction(Chlamy,'Glucose-1-phosphate adenylyltransferase 4 y',{'ATP [y]','G1P [y]','AGPase [y]','AGPase_ADP_ADPG_complex [y]'},[1 1 1 -1],false);

Chlamy=addReaction(Chlamy,'ADPGstore y',{'ADPG [y]'},-1,true);  

%% Diffusion stroma - pyrenoid
Chlamy=addReaction(Chlamy,'Diffusion RuBP',{'RuBP [s]','RuBP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion 3PGA',{'3PGA [s]','3PGA [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion GAP',{'GAP [s]','GAP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion DHAP',{'DHAP [s]','DHAP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion FBP',{'FBP [s]','FBP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion F6P',{'F6P [s]','F6P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion E4P',{'E4P [s]','E4P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion X5P',{'X5P [s]','X5P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion SBP',{'SBP [s]','SBP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion S7P',{'S7P [s]','S7P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion R5P',{'R5P [s]','R5P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion Ru5P',{'Ru5P [s]','Ru5P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion BPGA',{'BPGA [s]','BPGA [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion G1P',{'G1P [s]','G1P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion G6P',{'G6P [s]','G6P [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion ADPG',{'ADPG [s]','ADPG [y]'},[-1 1],true);

Chlamy=addReaction(Chlamy,'FNR',{'NADP [s]','NADPH [s]'},[-1 1],false); 
Chlamy=addReaction(Chlamy,'ATPase',{'ADP [s]','ATP [s]'},[-1 1],false); 

Chlamy=addReaction(Chlamy,'Diffusion ADP',{'ADP [s]','ADP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion ATP',{'ATP [s]','ATP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion NADP',{'NADP [s]','NADP [y]'},[-1 1],true);
Chlamy=addReaction(Chlamy,'Diffusion NADPH',{'NADPH [s]','NADPH [y]'},[-1 1],true);

Chlamy=addReaction(Chlamy,'Diffusion DHAP hs-c',{'DHAP [s]'},-1,true);
% Chlamy=addReaction(Chlamy,'Diffusion DHAP hs-c r',{'DHAP [s]'},1,false);
Chlamy=addReaction(Chlamy,'Diffusion BPGA hs-c',{'BPGA [s]'},-1,true);
% Chlamy=addReaction(Chlamy,'Diffusion BPGA hs-c r',{'BPGA [s]'},1,false);
Chlamy=addReaction(Chlamy,'Diffusion GAP hs-c',{'GAP [s]'},-1,true);
% Chlamy=addReaction(Chlamy,'Diffusion GAP hs-c r',{'GAP [s]'},1,false);
Chlamy=addReaction(Chlamy,'Diffusion 3PGA hs-c',{'3PGA [s]'},-1,true);
% Chlamy=addReaction(Chlamy,'Diffusion 3PGA hs-c r',{'3PGA [s]'},1,false);

Chlamy_splitted = split_rxns(Chlamy);
% 
save('Chlamy_model.mat','Chlamy_splitted')
% Chlamy_sbml=convertCobraToSBML(Chlamy,3,1,{'s','y'},{'stroma','pyrenoid'})
% OutputSBML(Chlamy_sbml, 'Chlamy.xml');
import os

tikxlist = [
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxTNoneyyTNone.jsonxcomp.tex",
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxT0.0yyT0.0.jsonxcomp.tex",
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxT0.1yyT0.1.jsonxcomp.tex",
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxTNoneyyTNone.jsonycomp.tex",
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxT0.0yyT0.0.jsonycomp.tex",
    "NStokesNwtnit4_time1.0_nu0.01_mesh25_Nts40_dt0.025NY4NU4alphau1e-07yxT0.1yyT0.1.jsonycomp.tex"]

outnamelist = ["yxTNoneyyTNonex.pdf",
               "yxT00yyT00x.pdf",
               "yxT01yyT01x.pdf",
               "yxTNoneyyTNoney.pdf",
               "yxT00yyT00y.pdf",
               "yxT01yyT01y.pdf"]

for indx in range(len(tikxlist)):
    os.system("cp " + tikxlist[indx] + " pgfpictoexp.tex")
    os.system("pdflatex exptikx.tex")
    os.system("cp exptikx.pdf " + outnamelist[indx])

// This file contains functions which record data.

/*------------------- Recording the progress of simulations ------------------*/

void RecordProgressSubroutine(int *cntMe, FILE *fOut) {
  int *cntAll, sumCnt;

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt);
	Printf("%5d ", sumCnt);
  }
  free(cntAll);
}

void RecordProgressSubroutine2(int *nPar1, int *nPar2, int nPar3, 
		FILE *fOut, int mode) {
  int n, *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  V2SET_ALL(sumCnt, 0);
  for(n = 0; n < 2; n++) {
	CONT(!((n == 0 && (int)(mode / 2) == 1) || (n == 1 && mode % 2 == 1)));
	MPI_Gather(((n == 0) ? nPar1 : nPar2), 1, MPI_INT, cntAll, 
			1, MPI_INT, 0, MPI_COMM_WORLD);
	CONT(rank != 0);
	sumCnt[n] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[n]);
	Printf("%5d/", sumCnt[n]);
  }
  if (rank == 0) {
	if (mode == 0) { 
		fprintf(fOut, "%d\t", nPar3);
		Printf("%5d ", nPar3);
	}
	else {
		fprintf(fOut, "%d\t%d\t", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
		Printf("%5d/%5d ", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
	}
  }
  free(cntAll);
}

void RecordProgressSubroutine3(int *cntMe, int *cntMe2, FILE *fOut) {
  int *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[0] = SumArrInt(cntAll, nCpu);
  }
  MPI_Gather(cntMe2, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[1] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[0] + sumCnt[1]);
	Printf("%5d ", sumCnt[0] + sumCnt[1]);
  }
  free(cntAll);
}

// Record the progress of the simulation.
void RecordProgress(void) {
  int n, nAcp, mode, nMot2[3], *pArr;
  double vol, area;
  time_t now;
  FILE *fOut;

  if (rank == 0) {
	fOut = fopen(GenFileName("Progress"), "a");
	fprintf(fOut, "%lld\t%g\t", currTimeStep, (double)currTimeStep * dtReal);
	Printf("%12d ", currTimeStep);
  }
  nAcp = nAbp - nMot;
  // Gather information from all subdomains first.
  if (nAct > 0) {
	mode = 0;
	if (actAss.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&actM.c, &actM.c, nAct, fOut, mode);
	if (actNuc.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ RecordProgressSubroutine(&nActFilaMe, fOut); }
	if (actNuc.gTgl != 0) { RecordProgressSubroutine(&actNuc.cntMe, fOut); }
	if (actAss.gTgl != 0) { RecordProgressSubroutine(&actAss.cntMe, fOut); }
	if (actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ RecordProgressSubroutine(&actDis.cntMe, fOut); }
	if (actBst.gTgl != 0) { RecordProgressSubroutine(&actBst.cntMe, fOut); }
	if (actSev.gTgl != 0) { RecordProgressSubroutine(&actSev.cntMe, fOut); }
	if (actAnn.gTgl != 0) { RecordProgressSubroutine(&actAnn.cntMe, fOut); }
	if (actCap.gTgl != 0) { RecordProgressSubroutine(&actCap.cntMe, fOut); }
	if (actUnc.gTgl != 0) { RecordProgressSubroutine(&actUnc.cntMe, fOut); }
  }
  if (nAcp > 0) { 
	mode = 0;
	if (acpInaUnb.gTgl != 0 || acpMoBind.gTgl != 0) { mode += 2; }
	if (acpUnb.gTgl != 0 || acpReb.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&nAcpMme, &nAcpInaMe, nAcp, fOut, mode); 
	if (acpInaUnb.gTgl != 0) 
	{ RecordProgressSubroutine(&acpInaUnb.cntMe, fOut); }
	if (acpMoBind.gTgl != 0) 
	{ RecordProgressSubroutine(&acpMoBind.cntMe, fOut); }
	if (acpUnb.gTgl != 0) { RecordProgressSubroutine(&acpUnb.cntMe, fOut); }
	if (acpReb.gTgl != 0) { RecordProgressSubroutine(&acpReb.cntMe, fOut); }
  }
  if (nMot > 0) { 
	mode = 0;
	if (motInaUnb.gTgl != 0 || motMoBind.gTgl != 0) { mode += 2; }
	if (motUnb.gTgl != 0 || motReb.gTgl != 0) { mode++; }
	if (motSA.gTgl == 0) { 
		RecordProgressSubroutine2(&nMotMme, &nMotInaMe, nMot, fOut, mode); 
		if (motInaUnb.gTgl != 0) 
		{ RecordProgressSubroutine(&motInaUnb.cntMe, fOut); }
		if (motMoBind.gTgl != 0) 
		{ RecordProgressSubroutine(&motMoBind.cntMe, fOut); }
		if (motUnb.gTgl != 0) { RecordProgressSubroutine(&motUnb.cntMe, fOut); }
		if (motReb.gTgl != 0) { RecordProgressSubroutine(&motReb.cntMe, fOut); }
	}
	else {
		nMot2[0] = 0;
		FOR_ABPME(n) {
			pArr = &P2A(abp.ch,n,0,nChAb);
			CONT(pArr[2] != 2);
			CONT(!(pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0));
			nMot2[0] += 2;
		}
		nMot2[1] = 2 * nMotMme - nMot2[0] + nMotInaMe;
		nMot2[2] = nMot * 2;
		RecordProgressSubroutine2(&nMot2[0], &nMot2[1], nMot2[2], fOut, mode); 
		if (motUnb.gTgl != 0) {
			RecordProgressSubroutine3(&motUnb.cntMe, &motInaUnb.cntMe, fOut); 
		}
		if (motReb.gTgl != 0) { 
			RecordProgressSubroutine3(&motReb.cntMe, &motMoBind.cntMe, fOut); 
		}
	}
	if (motWalk.gTgl != 0) { RecordProgressSubroutine(&motWalk.cntMe, fOut); }
    if (motSA.gTgl != 0) {
		RecordProgressSubroutine(&motSA.cntNucMe, fOut);
		RecordProgressSubroutine(&motSA.cntAssMe, fOut);
		if (motSA.to.gTgl != 0) {
			RecordProgressSubroutine(&motSA.cntTurnMe, fOut);
		}
    }
  }
  if (rank == 0) {
	now = time(NULL);
	fprintf(fOut, "%s", ctime(&now));
	fclose(fOut);
	Printf("\n");
  }

}

/*------------------- Recording the progress of simulations ------------------*/

/*------------------------- Recording parameter values -----------------------*/

void RecordInitParameterSubroutine(const char *str, int tgl1, int tgl2, 
		int tgl3, FILE *fOut) {
  fprintf(fOut, "%s = %s / %s / %s\n", str, 
		((tgl1 != 0 && tgl3 != 0) ? "on" : "off"),
		((tgl1 != 0 && tgl2 != 0) ? "on" : "off"), (tgl1 != 0 ? "on": "off"));
}

void RecordInitParameterSubroutine2(const char *str, 
		FuncCont *cont, FILE *fOut) {
  fprintf(fOut, "%s = %s\n", str, (cont->tgl != 0) ? "on" : "off");
  if (cont->tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", cont->prdR, cont->prd);
  }
}

// Record the parameter values which were initially loaded from "condition" 
// and "Config".
void RecordInitParameter(void) {
  int chk[NDIM][2], chk2, k, nAcp;
  double dimDomC[NDIM] , vol, ratio, ratio2, cAbp, rAbp;
  char direc[4] = "xyz", ch[80];
  FILE *fOut;  

  nAcp = nAbp - nMot;

  V3DIV(dimDomC, dimDom, nCell);
  vol = V3PROD(dimDom) * CUBE(L_SCALE_IN_M) * 1e3;
  cAbp = (double)nAbp * 2. / N_AVO / vol * 1e6;
  rAbp = (double)nAbp / ((double)nAct * nActPerSeg);
  fOut = fopen(GenFileName("Parameter"), "w");

  fprintf(fOut, "Duration of initial network formation = %g [s] "
		"(%lld timesteps)\n", netForm.durR, netForm.dur);
  fprintf(fOut, "Time step = %g [s] (%g)\n", dtReal, dt);
  fprintf(fOut, "Threshold of unstable force = %g [pN] (%g)\n",
		MAG_UNSTABLE_FORCE * 1e12, magUnstF);

  fprintf(fOut, "\n================== Number and concentration of actins and "
		"ABPs =================\n");
  fprintf(fOut, "Number of actin segments = %d (C_A = %g uM)\n", nAct, 
		(double)nAct * 2. * nActPerSeg / N_AVO / vol * 1e6);
  fprintf(fOut, "Length of actin segments = %g [nm]\n", L_SCALE_IN_NM);
  fprintf(fOut, "Number of ABPs = %d (C_ABP = %g uM, R_ABP = %g)\n", nAbp, 
		cAbp, rAbp);
  if (nAbp > 0) {
	ratio = (nAbp > 0) ? (double)nMot / (double)nAbp : 0.;
	ratio2 = (nAcp > 0) ? (double)nAbpDet[0] / (double)nAcp : 0.;
	fprintf(fOut, "(N_ACP^C = %d, C_ACP^C = %g uM, R_ACP^C = %g)\n", nAbpDet[0],
			cAbp * (1. - ratio) * ratio2, rAbp * (1. - ratio) * ratio2);
	fprintf(fOut, "(N_ACP^B = %d, C_ACP^B = %g uM, R_ACP^B = %g)\n", nAbpDet[1],
			cAbp * (1. - ratio) * (1. - ratio2), rAbp * (1. - ratio) 
			* (1. - ratio2));
	fprintf(fOut, "(N_M = %d, C_M = %g uM, R_M = %g)\n", nMot, 
			cAbp * ratio, rAbp * ratio);
  }

  fprintf(fOut, "\n===================== Conditions for domain and boundaries "
		"=====================\n");
  fprintf(fOut, "Domain width: x = %g [um] (%g), y = %g [um] (%g),"
		" z = %g [um] (%g)\n", L_S2UM(dimDom[0]), dimDom[0], 
		L_S2UM(dimDom[1]), dimDom[1], L_S2UM(dimDom[2]), dimDom[2]);
  fprintf(fOut, "Periodic boundary condition: x = %s, y = %s, z = %s\n", 
		(pbc[0] == 1 ? "on" : "off"), (pbc[1] == 1 ? "on" : "off"),
		(pbc[2] == 1 ? "on" : "off"));
  fprintf(fOut, "Geometry of a domain = ");
  if (dir2D > -1) {
	fprintf(fOut, "2-D with a normal direction in %c\n", direc[dir2D]);
  }
  else {
	fprintf(fOut, "3-D\n");
  }
  fprintf(fOut, "Movement of boundaries following stress-strain relation"
		" = %s\n", ((bndMv.gTgl != 0) ? "on" : "off"));
  if (bndMv.gTgl != 0) {
	if (bndMv.stfUnit == 0) {
		fprintf(fOut, "Young's modulus = %g on -x, %g on +x, %g on -y, %g on "
				"+y, %g on -z, %g on +z [Pa]\n", S2PA(bndMv.stf[0]), 
				S2PA(bndMv.stf[3]), S2PA(bndMv.stf[1]), S2PA(bndMv.stf[4]), 
				S2PA(bndMv.stf[2]), S2PA(bndMv.stf[5]));
	}
	else {
		fprintf(fOut, "Young's modulus = %g on -x, %g on +x, %g on -y, %g on "
				"+y, %g on -z, %g on +z [N/m]\n", KS_S2NPM(bndMv.stf[0]), 
				KS_S2NPM(bndMv.stf[3]), KS_S2NPM(bndMv.stf[1]), 
				KS_S2NPM(bndMv.stf[4]), KS_S2NPM(bndMv.stf[2]), 
				KS_S2NPM(bndMv.stf[5]));
	}
	fprintf(fOut, "Conservation of volume encapsulated by boundaries = %s\n",
			((bndVol.gTgl != 0) ? "on" : "off"));
	if (bndVol.gTgl != 0) {
		fprintf(fOut, "Stiffnesses of volume conservation = %g [N/m^5] (%g)\n",
				bndVol.stf * (KT_IN_J / pow(L_SCALE_IN_M, 6.)), bndVol.stf);
	}
  }
  if (recBndLoc.prdR > 0.) {
	fprintf(fOut, "Period of recording the boundary location = %g [s] "
			"(%d timesteps)\n", recBndLoc.prdR, recBndLoc.prd);
  }
  if (nAct > 0) {
	fprintf(fOut, "Unbinding of actin filaments from boundaries"
	 		" = %s\n", ((bndUnb.gTgl != 0) ? "on" : "off"));
	if (bndUnb.gTgl != 0) { 
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				K0_BND_UNB);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_BND_UNB);
		fprintf(fOut, "Factor that adjusts k0 = %g on -x, %g on "
				"+x, %g on -y, %g on +y, %g on -z, %g on +z\n", bndUnb.facK0[0],
				bndUnb.facK0[1], bndUnb.facK0[2], bndUnb.facK0[3], 
				bndUnb.facK0[4], bndUnb.facK0[5]);
		fprintf(fOut, "Factor that adjusts x = %g on -x, %g on "
				"+x, %g on -y, %g on +y, %g on -z, %g on +z\n", bndUnb.facX[0], 
				bndUnb.facX[1], bndUnb.facX[2], bndUnb.facX[3], bndUnb.facX[4], 
				bndUnb.facX[5]);
	}
	fprintf(fOut, "Binding of actin filaments on boundaries"
			" = %s\n", ((bndReb.gTgl != 0) ? "on" : "off"));
	if (bndReb.gTgl != 0) {
		fprintf(fOut, "Binding of entire parts of the filaments"
				" = %s\n", ((bndReb.gTglPa != 0) ? "on" : "off"));
		fprintf(fOut, "Thickness of space at which filaments can bind = %g "
				"[nm] (%g)\n", bndReb.depR, bndReb.dep);
		fprintf(fOut, "Binding rate(k) = %g [1/uM s]\n", K_BND_REB);
		fprintf(fOut, "Factor that adjusts k = %g on -x, %g on +x, %g on -y, "
				"%g on +y, %g on -z, %g on +z\n", bndReb.facK[0], 
				bndReb.facK[1], bndReb.facK[2], bndReb.facK[3], bndReb.facK[4],
				bndReb.facK[5]);
	}
	if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
		fprintf(fOut, "Manuration rate of links between actin filaments and "
				"boundary(k0) = %g [1/s]\n", K0_ACT_MAT);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_ACT_MAT);
		fprintf(fOut, "Factor that adjusts k0 = %g\n", bndMat.facK0);
		fprintf(fOut, "Factor that adjusts x = %g\n", bndMat.facX);
		fprintf(fOut, "Maximum level of maturation = %d\n", bndMat.maxNFA);
	}
	fprintf(fOut, "Moving of actin filaments on boundaries"
			" = %s\n", ((bnd.gTglActMv != 0) ? "on" : "off"));
	if (bnd.gTglActMv != 0) {
		fprintf(fOut, "Drag coefficient for the movement (relative to actin "
				"segment) = %g on -x, %g on +x, %g on -y, %g on +y, %g on -z, "
				"%g on +z\n", bnd.drag[0], bnd.drag[1], bnd.drag[2], 
				bnd.drag[3], bnd.drag[4], bnd.drag[5]);
	}
	if (bndUnb.gTgl != 0 || bndReb.gTgl != 0) {
		if (recBndUnbReb.prdR > 0.) {
			fprintf(fOut, "Period of recording the unbinding and binding on "
					"boundaries = %g [s] (%d timesteps)\n", recBndUnbReb.prdR, 
					recBndUnbReb.prd);
		}
	}
  }

  fprintf(fOut, "\n======================== Information about measurement "
			"=========================\n");
  fprintf(fOut, "Net duration of simulation = %g [s] "
		"(%lld timesteps)\n", rheo.durR, rheo.dur);
  fprintf(fOut, "Duration of motor activation before rheological measurement "
		"= %g [s] (%lld timesteps)\n", motActiv.durR, motActiv.dur);

  fprintf(fOut, "\n============= Conditions for dynamic behaviors of actins "
		"and ABPs ==============\n");
  RecordInitParameterSubroutine2("Equilibration of counter for monomeric "
		"actins and ABPs", &updMono, fOut);

  if (nAct > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors of "
			"actins ==================\n");
	RecordInitParameterSubroutine("Nucleation of actins", 
			actNuc.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	fprintf(fOut, "Control of the number of actin filaments by the "
			"nucleation = %s\n", (actNuc.gTglFN != 0 ? "on" : "off"));
	RecordInitParameterSubroutine("Assembly of actins", 
			actAss.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Disassembly of actins", 
			actDis.gTgl, gTglActDynPres, gTglActDynNF, fOut);
  }
  if (nAcp > 0) {
	fprintf(fOut, "\n=================== Conditions for dynamic behaviors "
			"of ACPs ===================\n");
	RecordInitParameterSubroutine("Unbinding of inactive ACPs", 
			acpInaUnb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of monomeric ACPs", 
			acpMoBind.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Unbinding of active ACPs", 
			acpUnb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of inactive ACPs", 
			acpReb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	fprintf(fOut, "Check cross-linking angle at binding of inactive ACPs = "
			"%s\n", (acpReb.gTglCrsAng != 0 ? "on" : "off"));
	fprintf(fOut, "Implicit consideration of monomeric ACPs = %s\n", 
			(gTglImpAcpM != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors "
			"of motors ==================\n");
	RecordInitParameterSubroutine("Walking of active motors", 
			motWalk.gTgl, gTglMotWalkPres, gTglMotWalkNF, fOut);
	fprintf(fOut, "Slide-off by the walking = %s\n",
			(gTglMotWalkSld != 0 ? "on" : "off"));
	fprintf(fOut, "Self-assembly of motors = %s\n",
			(motSA.gTgl != 0 ? "on" : "off"));
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Regulate relative constant size of the self-assembled "
				"structure = %s\n",	(motSA.gTglConSiz != 0 ? "on" : "off"));
		RecordInitParameterSubroutine("Unbinding of motors", 
				motUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of motors", 
				motReb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Check the alignment between actin and motor filaments "
				"for binding of motors = %s\n", (motReb.gTglCrsAng != 0 
				? "on" : "off"));
		fprintf(fOut, "Allow the binding of motors in a dimer only in opposite "
				"directions = %s\n", (motReb.gTglOppDir != 0 ? "on" : "off"));
		fprintf(fOut, "Turnover of motor filaments = ");
		if (motSA.to.gTgl == 0) {
			fprintf(fOut, "no\n");
		}
		else {
			if (motSA.to.mode == 0) {
				fprintf(fOut, "free\n");
			}
			else {
				fprintf(fOut, "free & force-dependent\n");
			}
		}
  		RecordInitParameterSubroutine2("Equilibration of counter for motor "
				"self-assembly", &updMotN, fOut);
	}
	else {
		RecordInitParameterSubroutine("Unbinding of inactive motors", 
				motInaUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of monomeric motors", 
				motMoBind.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Unbinding of active motors", 
				motUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of inactive motors", 
				motReb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Check cross-linking angle for binding of inactive "
				"motors = %s\n", (motReb.gTglCrsAng != 0 ? "on" : "off"));
		fprintf(fOut, "Implicit consideration of monomeric motors = %s\n", 
				(gTglImpMotM != 0 ? "on" : "off"));
	}
  }
  fprintf(fOut, "\n======== Adjustment for rates of the dynamic behaviors of "
		"actin and ABP ========\n");
  fprintf(fOut, "During network formation(i*k) = %g\n", netForm.facK);
  fprintf(fOut, "During prestrain/prestress(i*k) = %g\n", pres.facK);
  fprintf(fOut, "During motor activation(i*k) = %g\n", motActiv.facK);

  if (nAct > 0) {
	if (actNuc.gTgl != 0 || actAss.gTgl != 0 || actDis.gTgl != 0 
			|| actBst.gTgl != 0 || actSev.gTgl != 0 || actAnn.gTgl != 0 
			|| actCap.gTgl != 0 || actUnc.gTgl != 0 || actDgd.gTgl != 0) {
		fprintf(fOut, "\n======================== Parameters of actin dynamics "
				"=========================\n");
	}
	if (actNuc.gTgl != 0) {
		fprintf(fOut, "Nucleation rate = %g [1/uM s]\n", actNuc.k);
	}
	if (actAss.gTgl != 0) {
		fprintf(fOut, "Assembly rate = %g [1/uM s] at barbed, %g [1/uM s] at "
				"pointed\n", actAss.k[0], actAss.k[1]);
	}
  }

  if (nAcp > 0) {
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0 || acpInaUnb.gTgl != 0 
			|| acpUnb.gTgl != 0) {
		fprintf(fOut, "\n================= Parameters for unbinding and "
				"binding of ACPs =================\n");
	}
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0) {
		fprintf(fOut, "Binding and binding rate(k) = %g [1/uM s]\n", 
				K_ACP_BIND);
		fprintf(fOut, "Factor that adjusts k = %g\n", acpReb.facK);
	}
	if (acpInaUnb.gTgl != 0 | acpUnb.gTgl != 0) {
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				K0_ACP_UNB);
	}
	if (acpUnb.gTgl != 0) {
		fprintf(fOut, "Factor that adjusts k0 = %g\n", acpUnb.facK0);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_ACP_UNB);
		fprintf(fOut, "Factor that adjusts x = %g\n", acpUnb.facX);
	}
  }

  if (nMot > 0) {
	if (motSA.gTgl != 0) { 
		fprintf(fOut, "\n====== Parameters for motor self-assembly, "
				"unbinding, binding, and walking =====\n");
	}
	else {
		fprintf(fOut, "\n============== Parameters for motor unbinding, "
				"binding, and walking ============\n");
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0) {  
		fprintf(fOut, "Binding rate(k) = %g [1/uM s]\n", 
				40. * motMC.nHead);
		fprintf(fOut, "Factor that adjusts k = %g\n", motReb.facK);
	}
	if (motInaUnb.gTgl != 0 || motUnb.gTgl != 0) {
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				(motUnb.gTgl != 0 || motInaUnb.gTgl != 0) ? 
				log(1 - motUnb.p[motUnb.maxF[0] - 1]) / REVSIGN(dtReal) : 0.);
	}
	if (motWalk.gTgl != 0) {
		fprintf(fOut, "Walking rate without applied force(k0) = %g [1/s]\n", 
				(motWalk.gTgl != 0) ? log(1 - motWalk.p[motWalk.maxF[0] - 1]) 
				/ REVSIGN(dtReal) : 0.);
		fprintf(fOut, "Factor that adjusts k0 of unbinding and walking = %g\n", 
				motUnb.facK0);
		fprintf(fOut, "Each walking distance = %g [nm]\n", 
				L_SCALE_IN_NM / (double)nChAcX);
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0 || motInaUnb.gTgl != 0 
			|| motUnb.gTgl != 0 || motWalk.gTgl != 0) {
		fprintf(fOut, "Transition rates between mechanochemical states in each "
				"head(k01, k10, k12, k21, k20) = %g, %g, %g, %g, %g [1/s]\n", 
				motMC.k01, motMC.k10, motMC.k12, motMC.k21, motMC.k20);
		fprintf(fOut, "Number of heads which each motor arm represents = %d\n", 
				motMC.nHead);
	}
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Average number of motors per each self-assembled "
				"structure = %d\n", motSA.nMotPerTF);
		fprintf(fOut, "k for motor assembly = %g [1/s]\n", motSA.kAss);
		fprintf(fOut, "Turnover of free motor filaments = %s\n", 
				(motSA.to.gTgl != 0) ? "on" : "off");
		if (motSA.to.gTgl != 0) {
			fprintf(fOut, "k for motor free turnover = %g [1/s]\n", motSA.to.k);
			if (motSA.to.mode == 1) {
				fprintf(fOut, "k0 for motor force-dependent turnover = %g "
						"[1/s]\n", motSA.to.k0);
				fprintf(fOut, "Sensitivity of motor force-dependent turnover "
						"to force = %g [m]\n", motSA.to.x);
			}
		}
	}
  }

  if (nAbp > 0) {
	fprintf(fOut, "\n============================== Geometry of ABP"
			" =================================\n");
	fprintf(fOut, "Number of binding sites for ABPs on each actin segment in "
			"longitudinal and transverse directions = %d, %d\n", nChAcX, 
			nChAcY);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Length of ACP^C arm = %g [nm] (%g)\n", 
			L_ACPC_ARM * 1.0e9, L_M2S(L_ACPC_ARM));
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Length of ACP^B arm = %g [nm] (%g)\n", 
			L_ACPB_ARM * 1.0e9, L_M2S(L_ACPB_ARM));
  }
  if (nMot > 0) {
	fprintf(fOut, "Length of motor arm = %g [nm] (%g)\n", 
			L_MOT_ARM * 1.0e9, L_M2S(L_MOT_ARM));
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Length of bare zone in self-assembled motors = %g [nm] "
				"(%g)\n", L_MOTBACK_CEN_DIST * 1.0e9, 
				L_M2S(L_MOTBACK_CEN_DIST));
		fprintf(fOut, "Spacing between motors in self-assembled motors = %g "
				"[nm] (%g)\n", L_MOTBACK_DIST * 1.0e9, L_M2S(L_MOTBACK_DIST));
	}
  }
  fprintf(fOut, "\n=============================== Drag coefficients "
		"==============================\n");
  fprintf(fOut, "Viscosity of medium = %g [Pa s] (%g)\n", VISCOSITY, 
		VISCOSITY * L_SCALE_IN_M / actF.dragR);
  if (nAct > 0) {
	fprintf(fOut, "Drag coefficient of actin segment = %g [kg/s] (%g)\n", 
			actF.dragR, 1.0);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^C = %g [kg/s] (%g)\n", 
			abpF.drag[0].n * actF.dragR, abpF.drag[0].n);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^B = %g [kg/s] (%g)\n", 
			abpF.drag[1].n * actF.dragR, abpF.drag[1].n);
  }
  if (nMot > 0) {
	fprintf(fOut, "Drag coefficient of motor = %g [kg/s] (%g)\n", 
			abpF.drag[2].n * actF.dragR, abpF.drag[2].n);
  }

  fprintf(fOut, "\n============================== Thermal fluctuation "
		"=============================\n");
  if (nAct > 0) {
	fprintf(fOut, "Thermal fluctuation of actin filaments = %s\n",
			(gTglActTherm != 0 ? "on" : "off"));
  }
  if (nAcp > 0) {
	fprintf(fOut, "Thermal fluctuation of ACPs = %s\n", 
			(gTglAcpTherm != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "Thermal fluctuation of motors = %s\n", 
			(gTglMotTherm != 0 ? "on" : "off"));
  }

  if (nAct > 0 || nAbp > 0) {
	fprintf(fOut, "\n========================= Stiffness of force models "
			"============================\n");
  }
  if (nAbp > 0) {
  fprintf(fOut, "Strength of repulsive forces between ABPs and between ABP "
		"and actin = %g [N/m] (%g)\n",	KS_S2NPM(abpF.rep.stf), 
		abpF.rep.stf);
  }
  if (nAct > 0) {
	fprintf(fOut, "------------------------------------- actin "
			"------------------------------------\n");
	fprintf(fOut, "Strength of repulsive forces between actin filaments = %g "
			"[N/m] (%g)\n",	KS_S2NPM(actF.rep.stf), actF.rep.stf);
	fprintf(fOut, "Bending stiffness of actin filaments = %g [Nm] (%g)\n",
			KB_S2NM(actF.bend.stf), actF.bend.stf);
	fprintf(fOut, "Persistence length of actin filaments = %g [um]\n", 
			KB_S2NM(actF.bend.stf) * L_SCALE_IN_M / KT_IN_J * 1e6);
	fprintf(fOut, "Extensional stiffness of actin filaments = %g [N/m] (%g)\n",
			KS_S2NPM(actF.spr.stf), actF.spr.stf);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "------------------------------------- ACP^C "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^C (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[0].stf), abpF.bend[0].stf);
	fprintf(fOut, "Bending stiffness of ACP^C (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[0].stf), abpF.a90[0].stf);
	fprintf(fOut, "Extensional stiffness of ACP^C = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[0].stf), abpF.spr[0].stf);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "------------------------------------- ACP^B "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^B (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[1].stf), abpF.bend[1].stf);
	fprintf(fOut, "Bending stiffness of ACP^B (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[1].stf), abpF.a90[1].stf);
	fprintf(fOut, "Extensional stiffness of ACP^B = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[1].stf), abpF.spr[1].stf);
  }
  if (nMot > 0) {
	fprintf(fOut, "------------------------------------- motor "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of motor (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[2].stf), abpF.bend[2].stf);
	fprintf(fOut, "Bending stiffness of motor (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[2].stf), abpF.a90[2].stf);
	fprintf(fOut, "Extensional stiffness 1 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf), abpF.spr[2].stf);
	fprintf(fOut, "Extensional stiffness 2 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf2), abpF.spr[2].stf2);
	if (motSA.gTgl != 0) {
	  	fprintf(fOut, "Bending stiffness of motor backbone = %g [Nm] (%g)\n",
				KB_S2NM(motSA.bend.stf), motSA.bend.stf);
	  	fprintf(fOut, "Extensional stiffness of motor backbone = %g [N/m] "
				"(%g)\n", KS_S2NPM(motSA.spr.stf), motSA.spr.stf);
	}
  }
  fprintf(fOut, "\n============================ Parallel processing "
			"===============================\n");
  fprintf(fOut, "Number of cores used for the job = %d\n", nCpu);
  fprintf(fOut, "Number of subdomains: %d in x, %d in y, %d in z\n", nCell[0], 
		nCell[1], nCell[2]);
  fprintf(fOut, "Initial Subdomain size: x = %g [um] (%g), y = %g "
		"[um] (%g), z = %g [um] (%g)\n", L_S2UM(dimDomC[0]), 
		dimDomC[0], L_S2UM(dimDomC[1]), dimDomC[1], L_S2UM(dimDomC[2]), 
		dimDomC[2]);
  RecordInitParameterSubroutine2("Adjust the size of subdomains to make CPU "
		"loads even", &updSubdSize, fOut);
  fprintf(fOut, "Thickness of overlapping regions = %g [um] (%g)\n",
		L_S2UM(neiEdge), neiEdge);
  FOR_NDIM(k) {
	chk[k][0] = (nCell[k] > 1) ? 1 : 0; 
	chk[k][1] = (nCell[(k+1)%NDIM] > 1 && nCell[(k+2)%NDIM] > 1) ? 1 : 0;
  }
  chk2 = (nCell[0] > 1 && nCell[1] > 1 && nCell[2] > 1)  ? 1 : 0;
  fprintf(fOut, "The ratio of an overlapping region per to a subdomain = %g "
		"[percent]\n", (2. * neiEdge * (dimDomC[0] * dimDomC[1] * chk[2][0]
		+ dimDomC[1] * dimDomC[2] * chk[0][0] + dimDomC[2] * dimDomC[0] * 
		chk[1][0]) - 4. * SQR(neiEdge) * (dimDomC[0] * chk[0][1] + dimDomC[1] 
		* chk[1][1] + dimDomC[2] * chk[2][1]) + 8. * CUBE(neiEdge) * chk2) 
		/ V3PROD(dimDomC) * 100.);	

  fprintf(fOut, "\n================================= Data loading "
		"=================================\n");
  fprintf(fOut, "Load network data from Config = %s\n",
		(gTglLoadNetData == 1 ? "on": "off"));
  fprintf(fOut, "Load information of actins which are previously fixed and "
		"clamped = %s\n", (gTglLoadNetDataFix == 1 ? "on": "off"));

  fprintf(fOut, "\n=============================== Data recording "
		"=================================\n");
  fprintf(fOut, "Directory for saving data = %s\n", dataFold);
  fprintf(fOut, "Deletion of pre-existing data files = %s\n", 
		(DELETE_FILE != 0 ? "on" : "off"));
  RecordInitParameterSubroutine2("Record Output and Progress", &recProg, fOut);

  fprintf(fOut, "----------------------------- Structural information "
		"---------------------------\n");
  RecordInitParameterSubroutine2("Record Config", &recConf, fOut);
  fprintf(fOut, "Recording structural information for visualization via VMD "
		"= %s\n", (recConfVmd.tgl != 0) ? "on" : "off");
  if (recConfVmd.tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", recConfVmd.prdR, 
			recConfVmd.prd);
	if (recConfVmd.mode == 0) { fprintf(fOut, "- Multiple files\n"); }
	else { fprintf(fOut, "- Single file\n"); }
	fprintf(fOut, "- Show boundaries in the network drawn via VMD = %s\n",
			(recConfVmd.gTglBnd != 0 ? "on" : "off"));
	fprintf(fOut, "- Record information for coloring a network in VMD = %s\n",
			(recConfVmd.gTglInfo != 0 ? "on" : "off"));
  }
  RecordInitParameterSubroutine2("Record the length of actin filaments", 
		&recFilaL, fOut);
  RecordInitParameterSubroutine2("Record distances between active ABPs", 
		&recCrsDist, fOut);

  fprintf(fOut, "--------------------------- Force, stress, and energy "
		"--------------------------\n");
  RecordInitParameterSubroutine2("Record longitudinal forces acting on ABPs", 
		&recLongF, fOut);
  RecordInitParameterSubroutine2("Record the mechanical energy of a network", 
		&recE, fOut);
  RecordInitParameterSubroutine2("Record internal elastic/viscous stresses",
		&recSecStre, fOut);
  if (recSecStre.tgl != 0) {
	fprintf(fOut, "- The number of measurements for the internal stresses = %d "
			"in x, %d in y, %d in z\n", nSecStreDiv[0], nSecStreDiv[1], 
			nSecStreDiv[2]);
  }

  fprintf(fOut, "--------------------- Dynamic behaviors of actin and ABPs "
		"----------------------\n");
  RecordInitParameterSubroutine2("Record actin severing", &recActSev, fOut);
  RecordInitParameterSubroutine2("Record ABP unbinding", &recAbpUnb, fOut);
  RecordInitParameterSubroutine2("Record ABP binding", &recAbpBind, fOut);
  RecordInitParameterSubroutine2("Record ABP turnover", &recAbpTurn, fOut);
  RecordInitParameterSubroutine2("Record how many times each dynamic event "
		"of ABPs occurs", &recAbpDyn, fOut);

  fprintf(fOut, "---------------------------------- Miscellany "
		"----------------------------------\n");
  RecordInitParameterSubroutine2("Record the information in unit of "
		"filaments", &recInfo, fOut);
  fclose(fOut);
}

/*------------------------- Recording parameter values -----------------------*/

/*-------------------------- Recording network structure ---------------------*/

// Record network configuration. Most of the data are related to the positions
// and chain information of particles.
void RecordConfig(char *fn, int period) {
  int n, k, *nActAll, *nAbpAll, *nMbAll;
  double *rActOff, *rAbpOff, *rMbOff;
  ListInt parL;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(nMbAll,int,nCpu);
  MALLOC(rActOff,double,nActC*NDIM);
  MALLOC(rAbpOff,double,nAbpC*NDIM);
  MALLOC(rMbOff,double,nMbC*NDIM);

  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  MPI_Gather(&nMbMe, 1, MPI_INT, nMbAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  if (rank == 0) {
	if (SumArrInt(nActAll, nCpu) != nAct) {
		Printf("Error: the sum of collected nActMe is different from "
				"nAct!\n");
		exit(-1);
	}
	if (SumArrInt(nAbpAll, nCpu) != nAbp) {
		Printf("Error: the sum of collected nAbpMe is different from "
				"nAbp!\n");
		exit(-1);
	}
	
	fOut = fopen(fn, ((period > 0) ? "a" : "w")); 
	if (period > 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", (int)(currTimeStep / period));
	}
	fprintf(fOut, "nAct : %d\n", nAct);
	fprintf(fOut, "nAbp : %d\n", nAbp);
	fprintf(fOut, "nMb : %d\n", nMb);
	fprintf(fOut, "dimDom : %g, %g, %g\n", dimDom[0], 
			dimDom[1], dimDom[2]);
	fprintf(fOut, "nActPerSeg, nChAcX, nChAcY : %d, %d, %d\n\n", nActPerSeg, 
			nChAcX, nChAcY);
  }
  // Offset the positions of particles. This can be recovered later using the
  // information of rGrid below.
  FOR_ACTME(n) { 
	FOR_NDIM(k) {
		P2(rActOff,n,k) = P2(act.r,n,k) - rGrid[k][0];
	}
  }
  FOR_ABPME(n) { 
	FOR_NDIM(k) {
		P2(rAbpOff,n,k) = P2(abp.r,n,k) - rGrid[k][0];
	}
  }
  FOR_MBME(n) { 
	FOR_NDIM(k) {
		P2(rMbOff,n,k) = P2(memb.r,n,k) - rGrid[k][0];
	}
  }
  // Positions and chain information
  if (rank == 0) { fprintf(fOut, "## Position for actin ##\n"); }
  RecordGather2dArrayDouble(nActMe, nActAll, 3, act.id, rActOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Position for ABP ##\n"); }
  RecordGather2dArrayDouble(nAbpMe, nAbpAll, 3, abp.id, rAbpOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for actin ##\n"); }
  RecordGather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for ABP ##\n"); }
  RecordGather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Position for membrane ##\n"); }
  RecordGather2dArrayDouble(nMbMe, nMbAll, 3, memb.id, rMbOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for membrane ##\n"); }
  RecordGather2dArrayInt(nMbMe, nMbAll, nChMb, memb.id, memb.ch, fOut, 0);
  if (rank == 0) {
	fprintf(fOut, "currTimeStep = %lld\n", currTimeStep);
	fprintf(fOut, "rGrid = %g, %g, %g, %g, %g, %g\n", rGrid[0][0], rGrid[1][0],
			rGrid[2][0], rGrid[0][nGrid[0] - 1], rGrid[1][nGrid[1] - 1], 
			rGrid[2][nGrid[2] - 1]);
	fprintf(fOut, "pbc = %d, %d, %d\n", pbc[0], pbc[1], pbc[2]);
	fprintf(fOut, "rheoWay = %d\n", rheoWay);
	fprintf(fOut, "bulkRheoWay = %d\n", bulkRheoWay);
	fprintf(fOut, "bulkRheoType = %d\n\n", bulkRheoType);
  }
  if (rank == 0) { 
	fprintf(fOut, "stra.accDsp = 0\n\n");
	fprintf(fOut, "meaStrePar = 0\n\n");
	fprintf(fOut, "appStraPar = 0\n\n");
  }
  if (recTraj.tgl != 0) {
	  RecordGather1dArrayIntWoIndWoCnt(recTraj.actMe.c, recTraj.actMe.l, 
			fOut, "actTraj");
  }
  else { 
	if (rank == 0) {
		fprintf(fOut, "actTraj = 0\n\n");
	}
  }
  if (recTraj.tgl2 != 0) {
	RecordGather1dArrayIntWoIndWoCnt(recTraj.abpMe.c, recTraj.abpMe.l, 
			fOut, "abpTraj");
  }
  else { 
	if (rank == 0) {
		fprintf(fOut, "abpTraj = 0\n\n");
	}
  }

  if (rank == 0) { 
	MALLOC(fixAct,int,nAct);
	if (motSA.gTgl != 0) { MALLOC(abpMotId,int,nAbp); }
  }
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  if (motSA.gTgl != 0) {
	Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.mId, abpMotId);
  }
  if (rank == 0) {
	MALLOC(parL.l,int,nAct);
	parL.c = 0;
	FOR_ACT(n) {
		CONT(fixAct[n] < 0); 
		InsertElement1dArrayWoChk(parL.l, &parL.c, n);
	} 
	fprintf(fOut, "fixAct = %d\n", parL.c);
	for(n = 0; n < parL.c; n++) {
		fprintf(fOut, "%d\t%d\n", parL.l[n], fixAct[parL.l[n]]);
	} 
	free(parL.l);
	fprintf(fOut, "\n");

	parL.c = 0;
	if (motSA.gTgl != 0) {
		MALLOC(parL.l,int,nAbp);
		FOR_ABP(n) {
			CONT(abpMotId[n] < 0); 
			InsertElement1dArrayWoChk(parL.l, &parL.c, n);
		} 
	}
	fprintf(fOut, "abp.mId = %d\n", parL.c);
	for(n = 0; n < parL.c; n++) {
		fprintf(fOut, "%d\t%d\n", parL.l[n], abpMotId[parL.l[n]]);
	} 
	free(fixAct);
	if (motSA.gTgl != 0) { 
		free(parL.l);
		free(abpMotId); 
	}
	fclose(fOut); 
  }
  free(rActOff);
  free(rAbpOff);
  free(rMbOff);
  free(nActAll);
  free(nAbpAll);
  free(nMbAll);
}

int RecordConfigVmdSubroutine(double *r1, double *r2, int *ind, int *mul) {
  int k, adjInd, fac[NDIM];
  double dr[NDIM];

  V3COPY(fac, ind);
  V3SUB(dr, r1, r2);
  FOR_NDIM(k) {
	CONT(!(fabs(dr[k]) > dimDomH[k]));
	fac[k] += (dr[k] > dimDomH[k] ? 1 : -1);
	if (fac[k] >= mul[k]) { fac[k] = 0; }
	else if (fac[k] < 0) { fac[k] = mul[k] - 1; } 
  }
  V3IND_BACK_INT(adjInd, fac, mul);

  return adjInd;
}

// Record configuration for visualization via VMD. Positions and chain 
// information are recorded in the PDB and PSF formats. Actin, ACP, and motor
// have different SEGNAMEs for distinction. 
// mode = 0: bonds between actins and ABPs are drawn to the end points of 
//           actin segments. 
//        1: actin segments are divided to sub-segments, so the bonds are
//           drawn to actual points. This will cost more memory in VMD.
void RecordConfigVmd(int mode) {
  int m, n, k, ind, ind2[NDIM], ind3, ind4, ind5, cnt, curr, amp = 1;
  int nRec, nD, ch, *pArr, *pArr2, loc[2], nCh, indMb;
  int nActVmd, nAbpVmd, mul[NDIM];
  int facInfoConfVmd, facBndConfVmd;
  int *nActAll, *nAbpAll, *chActVmd, *chAbpVmd, *iFila, *chkBond;
  int nMbVmd, *nMbAll, *chMbVmd, *rebMb, *fixMb, *chActAlt;
  int nBeadVmd;
  int bndCh[] = {0, 4, 1, 5, 2, 6, 3, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 
				2, 1, 3, 4, 6, 5, 7};
  double temp, dr[NDIM], domCen[NDIM], rPos[NDIM], sft[NDIM];
  double *rActVmd, *rAbpVmd, *occu, *mass;
  double *recActG, *recAbpG, *sendRecAct, *sendRecAbp, *pntRec;
  double *rMbVmd, *rBeadVmd;
  char **segType, fn[80], fnIn[80];
  FILE *fOut;
  ListInt chain;

  if (recConfVmd.mode == 0) {
	sprintf(fnIn, "ConfVmd_%d", (int)(currTimeStep / recConfVmd.prd));
  }
  else {
	sprintf(fnIn, "ConfVmd");
  }

  facInfoConfVmd = (recConfVmd.gTglInfo > 0) ? 1 : 0;
  facBndConfVmd = (recConfVmd.gTglBnd != 0) ? 8 : 0;

  MALLOC(chActAlt,int,nChAc);
  // If a direction has periodic boundary condition, the configuration is
  // automatically duplicated in the direction for a clearer view.
  FOR_NDIM(k) {
	mul[k] = (pbc[k] == 1) ? 2 : 1;
  }
  V3SET_ALL(mul, 1);

  nActVmd = nAct * V3PROD(mul) * ((mode != 0) ? nChAcX : 1);
  nAbpVmd = nAbp * V3PROD(mul);
  nMbVmd = 0;
  nBeadVmd = 0;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);

  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(iFila,int,nAct);
  } 
  // Gather the positions and chain information of particles.
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);  

  // Record a separate data file (.dat) containing information for coloring 
  // methods via VMD.
  nD = 6;
  if (recConfVmd.gTglInfo > 0) {
	MALLOC(sendRecAct,double,nActC * nD);
	MALLOC(sendRecAbp,double,nAbpC * nD);
	if (rank == 0) {
		MALLOC(recActG,double,nAct * nD);
		MALLOC(recAbpG,double,nAbp * nD);
	}
	FOR_ACTME(n) { 
		ch = BinaryPackActinChainArray(act.id[n]);
		V6SET(&P2A(sendRecAct,n,0,nD), recAct.len[n], 
				P2A(recAct.allF,n,NDIM,NDIM + 1), recAct.sprF[n], 
				recAct.bendF[n], ch, (double)recAct.cnt[n]);
	}
	FOR_ABPME(n) { 
		ch = BinaryPackAbpChainArray(abp.id[n]);
		V6SET(&P2A(sendRecAbp,n,0,nD), P2A(recAbp.len,n,0,recAbp.nL), 
				P2A(recAbp.allF,n,NDIM,NDIM + 1), recAbp.sprF[n], 
				recAbp.bendF[n], ch, (double)recAbp.cnt[n]);
	}
	Gather2dArrayDouble(nActMe, nActAll, nD, act.id, sendRecAct, recActG);	
	Gather2dArrayDouble(nAbpMe, nAbpAll, nD, abp.id, sendRecAbp, recAbpG);	

	if (rank == 0) {
		for(m = 0; m < 2; m++) {
			nRec = (m == 0) ? nAct : nAbp;
			pntRec = (m == 0) ? recActG : recAbpG;
			for(n = 0; n < nRec; n++) {
				for(k = 0; k < 4; k++) {
					if (P2A(pntRec,n,nD - 1,nD) == 0) { 
						V4SET_ALL(&P2A(pntRec,n,0,nD), 0.);
					}
					else {
						P2A(pntRec,n,k,nD) /= P2A(pntRec,n,nD - 1,nD);
		 				if (k > 0) {
							P2A(pntRec,n,k,nD) = F_S2PN(P2A(pntRec,n,k,nD));
						}
					}
				}
			}
		}
	}
	free(sendRecAct);
	free(sendRecAbp);
  }

  // Prepare the (duplicated) positions and chain information of particles.
  if (rank == 0) {
	MALLOC(rActVmd, double, nActVmd * NDIM);
	MALLOC(rAbpVmd, double, nAbpVmd * NDIM);
	MALLOC(chActVmd, int, nActVmd * ((mode != 0) ? nChAcY + 2 : nChAc));
	MALLOC(chAbpVmd, int, nAbpVmd * nChAb);
	MALLOC(occu, double, nAct + nAbp);
	MALLOC(mass, double, nAct + nAbp);
	MALLOC2(segType, char, nActVmd + nAbp); 
	MALLOC(chain.l, int, (nActVmd + nAbpVmd * 4 + nMbVmd 
			* ((dimMbNuc == 2) ? 2 : 6)  
			+ ((recConfVmd.gTglBnd != 0) ? 12 : 0)) * 2);
	for (k = 0; k < nActVmd + nAbp; k++) {
		MALLOC(segType[k], char, 4);
	}
	memset(chActVmd, -1, sizeof(int) * nActVmd 
			* ((mode != 0) ? nChAcY + 2 : nChAc));
	memset(chAbpVmd, -1, sizeof(int) * nAbpVmd * nChAb);
	// Copy position and chain information to multiplied array
	for(ind = 0; ind < V3PROD(mul); ind++) {
		V3IND_ASSIGN_INT(ind, mul, 1, ind2);
		FOR_ACT(n) {
			pArr = &P2A(chAct,n,0,nChAc);
			memset(chActAlt, -1, sizeof(int) * nChAc);
			for (m = 0; m < nChAc; m++) {
				ind5 = pArr[m];
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAct,n,0), (m < 2 ? 
						&P2(rAct,ind5,0) : &P2(rAbp,ind5,0)), ind2, mul);
				chActAlt[m] = ind5 + ind4 * ((m < 2) ? nAct : nAbp);
			}
			if (mode == 0) {
				ind3 = n + ind * nAct;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				for (m = 0; m < nChAc; m++) {
					P2A(chActVmd,ind3,m,nChAc) = chActAlt[m];
				}
			}
			else {
				ind3 = (n + ind * nAct) * nChAcX;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				pArr2 = &P2A(chActVmd,ind3,0,nChAcY + 2);
				if (pArr[1] > -1) { 
					pArr2[1] = (chActAlt[1] + ind4 * nAct + 1) * nChAcX - 1;
				}
				if (pArr[0] > -1) {
					for(m = 1; m < nChAcX; m++) {
						// Position
						CalcPosOnActSeg(&P2(rAct,n,0), &P2(rAct,pArr[0],0),
								rPos, (double)m / (double)nChAcX);
						VSV3ADD(&P2(rActVmd,ind3 + m,0), rPos, dimDom, ind2);
						// Chain
						V2SET(&P2A(pArr2,m,0,nChAcY + 2), ind3 + m + 1, 
								ind3 + m - 1);
					}
					// Chain at ends
					pArr2[0] = ind3 + 1;
					P2A(pArr2,nChAcX - 1,0,nChAcY + 2) 
							= (chActAlt[0] + ind4 * nAct) * nChAcX;

					for(m = 2; m < nChAc; m++) {
						CONT(pArr[m] < 0);
						loc[0] = (m - 2) / nChAcY;
						loc[1] = (m - 2) % nChAcY;
						P2A(pArr2,loc[0],loc[1] + 2,nChAcY + 2) 
								= chActAlt[m] + ind4 * nAbp;
					}
				}
				else {
					for(m = 1; m < nChAcX; m++) {
						V3SET_ALL(&P2(rActVmd,ind3 + m,0), 0.);
					}
				}

			}
		}
		FOR_ABP(n) {
			ind3 = n + ind * nAbp;
			VSV3ADD(&P2(rAbpVmd,ind3,0),&P2(rAbp,n,0),dimDom,ind2);	
			for (m = 0; m < 2; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAct,ind5,0), ind2, mul);
				if (mode == 0) { 
					P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAct;
				}
				else {
					loc[0] = FindAbpActinChain(ind5, n, 1);
					loc[0] = (loc[0] - 2) / nChAcY;
					P2A(chAbpVmd,ind3,m,nChAb) 
							= (ind4 * nAct + ind5) * nChAcX + loc[0];
				}
			}
			P2A(chAbpVmd,ind3,2,nChAb) = P2A(chAbp,n,2,nChAb);
			for (m = 3; m < 4; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAbp,ind5,0), ind2, mul);
				P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAbp;
			}
		}
		FOR_MB(n) {
			ind3 = n + ind * nMb;
			VSV3ADD(&P2(rMbVmd,ind3,0),&P2(rMb,n,0),dimDom,ind2);	
			for (m = 0; m < ((dimMbNuc == 2) ? 2 : 6); m++) {
				CONT(!(P2A(chMb,n,m,nChMb) > -1));
				ind4 = RecordConfigVmdSubroutine(&P2(rMb,n,0), 
						&P2(rMb,P2A(chMb,n,m,nChMb),0), ind2, mul);
				P2A(chMbVmd,ind3,m,nChMb) = 
						P2A(chMb,n,m,nChMb) + ind4 * nMb;
			}
		}
		for(n = 0; n < bead.n; n++) {
			ind3 = n + ind * bead.n;
			VSV3ADD(&P2(rBeadVmd,ind3,0),&P2(bead.r,n,0),dimDom,ind2);	
		}
	}
	// Fill in segment type of actin, ACP, and motor.
	for(n = 0; n < nActVmd; n++) { 
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		pArr = &P2A(chAct,ind,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) {
			sprintf(segType[n], "G");
		}
		else if (pArr[0] > -1 && pArr[1] < 0) {
			sprintf(segType[n], "FP");
			if (mode == 1) {
				if (P2A(chActVmd,n,1,nChAcY + 2) > -1) {
					sprintf(segType[n], "F");
				}
			}
		}
		else if (pArr[0] < 0 && pArr[1] > -1) {
			sprintf(segType[n], "FB");
			if (mode == 1) {
				if (n % nChAcX > 0) {
					sprintf(segType[n], "G");
				}
			}
		}
		else { sprintf(segType[n], "F"); }
	}
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		cnt = 1;
		if (pArr[0] < 0) { cnt++; }
		if (pArr[1] < 0) { cnt++; }
		if (motSA.gTgl != 0) {
			if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) {
				cnt = 4;
			}
		}
		switch(pArr[2]) {
			case 0: sprintf(segType[n + nActVmd], "CR%d", cnt);	break;
			case 1: sprintf(segType[n + nActVmd], "BU%d", cnt);	break;
			case 2: sprintf(segType[n + nActVmd], "MO%d", cnt);	break;
		}
	}
	
	if (recConfVmd.gTglInfo > 0 && currTimeStep > 1) {
		for(m = 0; m < 2; m++) {
			nRec = (m == 0) ? nAct : nAbp;
			pntRec = (m == 0) ? recActG : recAbpG;
			for(n = 0; n < nRec; n++) {
				temp = TrimDblVal(P2A(pntRec,n,2,nD), recConfVmd.minF, 
						recConfVmd.maxF);
				mass[n + ((m == 1) ? nAct : 0)] 
						= 1. + 20. * (temp  - recConfVmd.minF) 
						/ (recConfVmd.maxF - recConfVmd.minF);
			}
		}
	}
	else {
		for(n = 0; n < nAct; n++) {
			mass[n] = 12.;
		}
		for(n = 0; n < nAbp; n++) {
			mass[n + nAct] = 14.;
		}
	}
	
	for(n = 0; n < nAct + nAbp; n++) { occu[n] = 1.; }
	if (recPerc.tgl != 0) {
		for(n = 0; n < nAct + nAbp; n++) {
			if (iPerc[n] / 10 == 1) { occu[n] = 2.; }
			if (iPerc[n] % 10 == 1) { occu[n] = 3.; }
		}
	}

	// Write PDB file
	sprintf(fn, "%s.pdb", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf (fOut, "REMARK\n");

	FOR_NDIM(k) {	
  		if (bndMv.gTgl != 0) {
			domCen[k] = 0.5 * ((2 - mul[k]) * P2A(rGridInit,0,k,NDIM)
					+ mul[k] * P2A(rGridInit,1,k,NDIM));
		}
		else {
			domCen[k] = 0.5 * mul[k] * dimDom[k];
		}
	}
	// actin
	for(n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
	    VV3SUB(&P2(rActVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  C   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + 1, P2(rActVmd,n,0) * amp, 
				P2(rActVmd,n,1) * amp, P2(rActVmd,n,2) * amp, 
				occu[ind], segType[n]);
	}
	// ABP
	for(n = 0; n < nAbpVmd; n++) {
		VV3SUB(&P2(rAbpVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  N   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + nActVmd + 1, P2(rAbpVmd,n,0) * amp, 
				P2(rAbpVmd,n,1) * amp, P2(rAbpVmd,n,2) * amp, 
				occu[(n % nAbp) + nAct], segType[(n % nAbp) + nActVmd]);
	}
	// If a coloring method is used, a reference particle is necessary.
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "ATOM %6d  O   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", nActVmd + nAbpVmd + 1, 0., 0., 0., 1., "REF");
	}
	// Visualization of boundaries using cylinders
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			V3IND_ASSIGN_CONST_INT(n, 2, ind2);
			V3SET_ALL(sft, 0.);
			FOR_NDIM(k) {
				rPos[k] = rGrid[k][(ind2[k] == 0) ? 0 : nGrid[k] - 1] 
						- domCen[k] + sft[k] 
						+ ((ind2[k] == 1) ? dimDom[k] * (mul[k] - 1) : 0.); 
			}
			fprintf(fOut, "ATOM %6d  H   LYS     1    %8.3f%8.3f%8.3f  %4.2f"
					"  0.00      %s\n", nAbpVmd + nActVmd + n + facInfoConfVmd 
					+ 1, rPos[0], rPos[1], rPos[2], 1., "BND");
		}
	}
	for(n = 0; n < nMbVmd; n++) {
	    VV3SUB(&P2(rMbVmd,n,0), domCen);
		for(k = 0; k < indMb; k++) {
			fprintf(fOut, "ATOM %6d  P   LYS     1    ", nAbpVmd + nActVmd 
					+ facInfoConfVmd + 1 + facBndConfVmd + n * indMb + k);
			for(m = 0; m < NDIM; m++) {
				if (dimMbNuc == 2) { 
					fprintf(fOut, "%8.3f", ((m == dirNormMbNuc) ? domCen[m] 
							* amp * (k * 2 - 1) :  P2(rMbVmd,n,m) * amp));
				}
				else {
					fprintf(fOut, "%8.3f", P2(rMbVmd,n,m) * amp);
				}
			}
			if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
				temp = (fixMb[n % nMb] > -1)  ? 2. : 1.;
			}
			else { temp = 1.; }
			fprintf(fOut, "  %4.2f  0.00      %s\n", temp, 
					(ISNUC(idxMb[n % nMb]) ? "NUC" : "MB")); 
		}
	}
	for(n = 0; n < nBeadVmd; n++) {
	    VV3SUB(&P2(rBeadVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  S   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", nActVmd + nAbpVmd + facInfoConfVmd + 1 
				+ facBndConfVmd + nMbVmd * indMb + n, P2(rBeadVmd,n,0) * amp, 
				P2(rBeadVmd,n,1) * amp, P2(rBeadVmd,n,2) * amp, 1., "BD");
	}

	fprintf(fOut, "ENDMDL\n");
	fclose(fOut);

	// Write PSF file
	sprintf(fn, "%s.psf", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf(fOut, "PSF CMAP\n\n");
	fprintf(fOut, "%8d !NTITLE\n\n", 7);
	fprintf(fOut, "%8d !NATOM\n", nActVmd + nAbpVmd + nMbVmd 
			* indMb + facInfoConfVmd + facBndConfVmd + nBeadVmd);
	// actin
	for (n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		fprintf(fOut, "%8d %-4s 1    LYS  C      11%11.5f%14.4f%11d\n", 
				n + 1, segType[n], ((iSupp[ind] > -1) ? 1. : 0.), 
				mass[ind], 0);
	}
	// ABP
	for (n = 0; n < nAbpVmd; n++) {
		fprintf(fOut, "%8d %-4s 1    LYS  N      38%11.5f%14.4f%11d\n", 
				n + nActVmd + 1, segType[(n % nAbp) + nActVmd], 
				((iSupp[(n % nAbp) + nAct] > -1) ? 1. : 0.),
				mass[(n % nAbp) + nAct], 0);
	}
	// A reference particle for coloring method
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "%8d %-4s 1    LYS  O      38%11.5f%14.4f%11d\n", 
				nAbpVmd + nActVmd + 1, "REF", 0., 16., 0);
	}
	// Visualization of boundaries
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			fprintf(fOut, "%8d %-4s 1    LYS  H      38%11.5f%14.4f%11d\n", 
					nAbpVmd + nActVmd + n + 1 + facInfoConfVmd, 
					"BND", 0., 1., 0);
		}
	}
	// Membrane
	for (n = 0; n < nMbVmd * indMb; n++) {
		if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
			cnt = (rebMb[(n / indMb) % nMb] > -1) ? 2 : 1;
			if (cnt == 2) {
				for(k = nChMb - nMbAct; k < nChMb; k++) {
					if (P2A(chMb,(n / indMb) % nMb,k,nChMb) > -1) { cnt++; }
				}
			}
		}
		else { cnt = 1; }
		fprintf(fOut, "%8d %-4s 1    LYS  P      38%11.5f%14.4f%11d\n", 
				nAbpVmd + nActVmd + 1 + facInfoConfVmd + facBndConfVmd + n,
				(ISNUC(idxMb[(n / indMb) % nMb]) ? "NUC" : "MB"), 
				0., (double)cnt, 0); 
	}
	for(n = 0; n < nBeadVmd; n++) {
	    VV3SUB(&P2(rBeadVmd,n,0), domCen);
		fprintf(fOut, "%8d %-4s 1    LYS  S      38%11.5f%14.4f%11d\n", 
				nActVmd + nAbpVmd + facInfoConfVmd + 1 + facBndConfVmd 
				+ nMbVmd * indMb + n, "BD", 0., 18., 0);
	}
	chain.c = 0;
	// Chain information from actin (chActVmd)
	nCh = (mode != 0) ? nChAcY + 2 : nChAc;
	for (n = 0; n < nActVmd; n++) {
		CONT(!(P2A(chActVmd,n,0,nCh) > -1 && P2A(chActVmd,n,1,nCh) < 0));
		curr = n;
		while(P2A(chActVmd,curr,0,nCh) > -1) {
			V3SUB(dr, &P2(rActVmd,curr,0), 
					&P2(rActVmd,P2A(chActVmd,curr,0,nCh),0));
			FOR_NDIM(k) {
				BREAK(pbc[k] != 0 && fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = curr + 1;
				P2A(chain.l,chain.c,1,2) = 
						P2A(chActVmd,curr,0,nCh) + 1;
				chain.c++; 
			}
			curr = P2A(chActVmd,curr,0,nCh);	
		}
	}
	// Chain information from ABP (chAbpVmd)
	for (n = 0; n < nAbpVmd; n++) {
		pArr = &P2A(chAbpVmd,n,0,nChAb);
		for(m = 0; m < nChAb - 1; m++) {
			CONT(m == 2);
			CONT(pArr[m] < 0);
			CONT(!(ISMTF(pArr[2])) && m > 2);
			if (m < 2) { 
				V3SUB(dr, &P2(rAbpVmd,n,0), &P2(rActVmd,pArr[m],0));
			}
			else {
				V3SUB(dr, &P2(rAbpVmd,n,0),	&P2(rAbpVmd,pArr[m],0));
			}
			FOR_NDIM(k) {
				BREAK(fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = n + nActVmd + 1;
				P2A(chain.l,chain.c,1,2) = pArr[m] + 1 
						+ ((m < 2) ? 0 : nActVmd);
				(chain.c)++; 
			}
		}
    }
	// Chain information is also required for cylinders showing boundaries.
	if (recConfVmd.gTglBnd != 0) {
		ind = nActVmd + nAbpVmd + 1 + facInfoConfVmd;
		for(n = 0; n < 12; n++) {
			P2A(chain.l,chain.c,0,2) = P2A(bndCh,n,0,2) + ind;
			P2A(chain.l,chain.c,1,2) = P2A(bndCh,n,1,2) + ind;
			(chain.c)++; 
		}	
	}
	// Chains for membrane
	ind = nActVmd + nAbpVmd + facInfoConfVmd + 1 + facBndConfVmd;
	for (n = 0; n < nMbVmd; n++) {
		for(m = 0; m < ((dimMbNuc == 2) ? 1 : 6); m++) {
			ind3 = P2A(chMbVmd,n,m,nChMb);
			CONT(!(ind3 > -1));
			if (dimMbNuc == 3) {
				CONT(P2A(chkBond,n,m,6) != 0);
			}
			V3SUB(dr, &P2(rMbVmd,n,0), &P2(rMbVmd,ind3,0));
			FOR_NDIM(k) {
				BREAK(fabs(dr[k]) > dimDomH[k])
			}
			CONT(k != NDIM);
			for(k = 0; k < indMb; k++) {
				P2A(chain.l,chain.c,0,2) = n * indMb + k + ind;
				P2A(chain.l,chain.c,1,2) = ind3 * indMb + k + ind;
				(chain.c)++; 
			}
			if (dimMbNuc == 3) {
		        loc[0] = FindElementArray(&P2A(chMbVmd,ind3,0,
						nChMb), 6, n, 0, 1);
				P2A(chkBond,ind3,loc[0],6) = 1;
			}
		}
	}
	fprintf(fOut, "\n%8d !NBOND: bonds\n", chain.c);
	for(n = 0; n < chain.c; n++) {
		fprintf(fOut, "%8d%8d", P2A(chain.l,n,0,2), P2A(chain.l,n,1,2));
		CONT((n + 1) % 4 != 0);
		fprintf(fOut, "\n"); 
	}
	if (chain.c % 4 != 0) {
		fprintf(fOut, "\n"); 
	}	
   	fclose(fOut);

	free(rActVmd);
	free(rAbpVmd);
	free(chActVmd);
	free(chAbpVmd);
	free(chain.l);
	free(occu);
	free(mass);
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
	free(iFila);
	for (n = 0; n < nActVmd + nAbp; n++) { free(segType[n]); }
	free(segType);
	if (recConfVmd.gTglInfo > 0) {
		free(recActG);
		free(recAbpG);
	}
  }
  free(nActAll);
  free(nAbpAll);
  free(chActAlt);
}

void RecordActinInfo(int period) {
  int m, n, k, ind, cnt, *pArr, *nActAll, *nAbpAll, *cntAct, *iFila;
  double tensF, *lAct;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(lAct,double,nAct);
	MALLOC(cntAct,int,nAct);
	MALLOC(iFila,int,nAct);
  }
  // Gather information about the positions and chains of actins 
  GatherActChainPosition(nActAll, chAct, rAct);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.len, lAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, recAct.cnt, cntAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);

  if (rank == 0) {
	fOut = fopen(GenFileName("ActInfo"), "a");
	if (currTimeStep == period) {
		fprintf(fOut, "%d\t0\t0\t0\t0\t0\t0\t0\n", nAct);
	}
	FOR_ACT(n) {
		tensF = (cntAct[n] == 0) ? 0. : F_S2PN(actF.spr.stf 
				* (lAct[n] / (double)cntAct[n] - actF.spr.eq));
		fprintf(fOut, "%d\t%g\t%g\t%g\t", n, P2(rAct,n,0), P2(rAct,n,1), 
				P2(rAct,n,2));
		fprintf(fOut, "%d\t%d\t%g\t%d\n", P2A(chAct,n,0,nChAc), 
				P2A(chAct,n,1,nChAc), tensF, iFila[n]);
	}
	fclose(fOut);
  }

  free(nActAll);
  if (rank == 0) {
	free(rAct);
	free(chAct);
	free(lAct);
	free(cntAct);
	free(iFila);
  }
}

void RecordConfigMatlab(void) {
}

/*-------------------------- Recording network structure ---------------------*/

/*--------------------- Structural properties of network ---------------------*/

// Record the distribution of actin filaments
void RecordFilamentLength(int period) {
  int n, len, curr, cntFila, maxLen, *filaLen, *nActAll, *pArr;
  double avgFilaLen, stdFilaLen, len2;
  FILE *fOut;

  maxLen = -1;
  MALLOC(nActAll,int,nCpu);
  if (rank == 0) { MALLOC(chAct,int,nAct*nChAc); }
  // Gather information about chains.
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct); 

  if (rank == 0) {
	MALLOC(filaLen,int,nAct+1);
	memset(filaLen, 0, sizeof(int) * (nAct + 1) );
	cntFila = 0; 
	avgFilaLen = 0.; 
	stdFilaLen = 0.;

	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) { 
			filaLen[0]++; 
		} 
		else if (pArr[0] > -1 && pArr[1] < 0) { 
			len = 1; 
			curr = pArr[0];
			while (P2A(chAct,curr,0,nChAc) > -1) {
				curr = P2A(chAct,curr,0,nChAc);
				len++;
			}
			filaLen[len]++;
			if (maxLen < len) { maxLen = len; }
		}
	}
	fOut = fopen(GenFileName("ActFilaLen"), "a");
	// Record the number of lines recorded.
    Fprintf1dArrayIntWFil(fOut, &maxLen, 1, 0, 3);
	// Calculate average length and standard deviation.
	for(n = 1; n < maxLen + 1; n++) {
		len2 = L_S2UM(n);
	    fprintf(fOut, "%d\t%g\t%d\n", n - 1, len2, filaLen[n]);
	    cntFila += filaLen[n];
		avgFilaLen += (double)(len2 * filaLen[n]);
		stdFilaLen += (double)(SQR(len2) * filaLen[n]);
	}
	if (cntFila == 0) { 
		avgFilaLen = 0.; 
		stdFilaLen = 0.;
	}
	else { 
		avgFilaLen /= (double)cntFila; 
		stdFilaLen = sqrt(stdFilaLen / (double)cntFila - SQR(avgFilaLen));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxLen - 1, avgFilaLen, stdFilaLen);
	fclose(fOut);
	free(filaLen);
	free(chAct);
  }
  free(nActAll);
}

// Record the distribution of crosslinking distance which corresponds to 
// distance between active ABPs on a single filament.
void RecordCrosslinkDistance(void) {
  int n, k, l, curr, len, maxLen, actInd, abpInd, cntDist, CS, cnt, prevIMot;
  int *nActAll, *nAbpAll, *crossDist, *pArr, *pArr2, *iMot, *iFila;
  double avgCrossDist, stdCrossDist, len2;
  FILE *fOut;
  ListInt motIFila;

  maxLen = NEG_LARGE_VALUE;
  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	if (motSA.gTgl != 0) {
		MALLOC(motIFila.l,int,nAbp*2);
		MALLOC(iFila, int, nAct);
		MALLOC(iMot, int, nAbp);
		memset(iMot, -1, sizeof(int) * nAbp);
	}
  }

  // Gather information about chains of actin
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);
  if (motSA.gTgl != 0) {
	Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  }

  if (rank == 0) {
	MALLOC(crossDist, int, nAct * nChAcX);
	fOut = fopen(GenFileName("CrossDist"),"a");
	memset(crossDist, 0, sizeof(int) * nAct * nChAcX);
	avgCrossDist = 0.;
	cntDist = 0;

	cnt = 0;
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		CONT(!(ISMTF(pArr[2])));
		CONT(!(pArr[3] > -1 && pArr[4] < 0));
		motIFila.c = 0;
		curr = n;
		while(curr > -1) { 
			for(k = 0; k < 2; k++) {
				actInd = P2A(chAbp,curr,k,nChAb);
				CONT(!(actInd > -1));
				InsertElement1dArrayWChk(motIFila.l, &motIFila.c, 
						iFila[actInd]);
			}
			curr = P2A(chAbp,curr,3,nChAb);
		}
		if (motIFila.c > 1) {
			curr = n;
			while(curr > -1) { 
				iMot[curr] = cnt;
				curr = P2A(chAbp,curr,3,nChAb);
			}	
			cnt++;
		}
	}
	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		CONT(!(pArr[0] > -1 && pArr[1] < 0));
		len = 0;
		CS = 0;
		prevIMot = -1;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = P2A(chAct,curr,2 + k * nChAcY + l,nChAc);
					CONT(abpInd < 0);
					pArr2 = &P2A(chAbp,abpInd,0,nChAb);
					if (ISMTF(pArr2[2])) {
						CONT(iMot[abpInd] < 0);
						if (iMot[abpInd] == prevIMot) {
							len = 0;
							continue;
						}
						prevIMot = iMot[abpInd];
					}
					else {
						CONT(!(pArr2[0] > -1 && pArr2[1] > -1));
						prevIMot = -1;
					}
					// if it is not the end of filament
					if (CS != 0) {
						crossDist[len]++;
						if (maxLen < len) { maxLen = len; }
					}
					len = 0;
					CS = 1;
				}
				len++;
			}
			curr = P2A(chAct,curr,0,nChAc);
		}
	}
	if (maxLen < 0) { maxLen = -1; }
	for (n = 0; n < maxLen + 1; n++) {
		cntDist += crossDist[n];
	}
	fprintf(fOut, "%d\t%d\t0\n", maxLen + 2, cntDist);
	// Calculate average and standard deviation.
	for (n = 0; n < maxLen + 1; n++) {
		len2 = L_S2NM((double)n / (double)nChAcX);
		fprintf(fOut, "%d\t%g\t%d\n", n, len2, crossDist[n]);
		avgCrossDist += (double)(crossDist[n] * len2);
		stdCrossDist += (double)(crossDist[n] * SQR(len2));
	}
	if (cntDist == 0) {
		avgCrossDist = 0.;
		stdCrossDist = 0.;
	}
	else { 
		avgCrossDist /= (double)cntDist; 
		stdCrossDist = sqrt(stdCrossDist / (double)cntDist - SQR(avgCrossDist));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxLen + 1, avgCrossDist, stdCrossDist);
	fclose(fOut);
	free(crossDist);
	free(chAct);
	free(chAbp);
	if (motSA.gTgl != 0) { 
		free(motIFila.l);
		free(iFila);
		free(iMot); 
	}
  }
  free(nActAll);
  free(nAbpAll);
}

/*--------------------- Structural properties of network ---------------------*/

/*----------------------------- Energy calculation ---------------------------*/

void RecordMechEnergy(int mode) {
  int m, n, k, locActInd, abpInd, *pArr, side, kind;
  double len, lenEq, ang, fac;
  double E[4], Esum[4], *Eall, dr[2][NDIM + 1];
  double rPos[NDIM], rPos2[NDIM], *rPnt[2];
  FILE *fOut;
  ListInt *pAbpL;
 
  if (rank == 0) { MALLOC(Eall,double,nCpu*4); }
  V4SET_ALL(E, 0.);

  // Bending of actin filament
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
    CONT(!(pArr[0] > -1 && pArr[1] > -1));
	CalcVec(dr[0], &P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0));
	CalcVec(dr[1], &P2(act.r,iAct[pArr[1]],0), &P2(act.r,n,0));
	ang = V3ANG(dr[0], dr[1]);
	E[0] += 0.5 * actF.bend.stf * SQR(ang);
  }

  // Extension of actin filament (double-counting)
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
	CONT(pArr[0] < 0);
	len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0), 0);
	E[1] += 0.5 * actF.spr.stf * SQR(len - actF.spr.eq);
  }

  // Bending of ABP (actin-ABP-actin)
  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
	kind = pArr[2];
	CONT(ISMTF(kind));
	for(k = 0; k < 2; k++) {
		dr[k][NDIM] = CalcVecDistActinAbp(dr[k], pArr[k], abp.id[n], 0);
	}
	V3REVSIGN(dr[1]);
	ang = V3ANG(dr[0], dr[1]);
	E[2] += 0.5 * abpF.bend[kind].stf * SQR(ang - abpF.bend[kind].eq);
  }
  // Bending of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		abpInd = abp.id[n];
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
	    CONT(!(pArr[3] > -1 && pArr[4] > -1));
		CalcVec(dr[0], &P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0));
		CalcVec(dr[1], &P2(abp.r,iAbp[pArr[4]],0), &P2(abp.r,n,0));
		ang = V3ANG(dr[0], dr[1]);
		E[2] += 0.5 * motSA.bend.stf * SQR(ang);
	}
  }

  // Extension of ABP
  FOR_ABPME(n){
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	kind = pArr[2];
	for(k = 0; k < 2; k++) {
		CONT(pArr[k] < 0);
		locActInd = iAct[pArr[k]];
		side = FindAbpActinChain(locActInd, abp.id[n], 0);
		rPnt[0] = &P2(act.r,locActInd,0);
		rPnt[1] = &P2(act.r,iAct[P2A(act.ch,locActInd,0,nChAc)],0);
		CalcPosOnActSegSide(rPnt[0], rPnt[1], rPos, side);
  		CalcVec(dr[0], rPos, &P2(abp.r,n,0));
		if (ISMTF(kind)) {
			CalcVec(dr[1], rPnt[0], rPnt[1]);
			fac = V3DOT(dr[0], dr[1]) / V3LEN_SQ(dr[1]);
			VS3SUB(rPos2, rPos, dr[1], fac);
			ApplyBoundCondVector(rPos2, -1, 0);
			CalcVec(dr[0], rPos2, &P2(abp.r,n,0));
			CalcVec(dr[1], rPos, rPos2);			
		}
		len = V3LEN(dr[0]);
		E[3] += 0.5 * abpF.spr[kind].stf 
				* SQR(len - abpF.spr[kind].eq);
		if (ISMTF(kind)) {	
			len = V3LEN(dr[1]);
			E[3] += 0.5 * abpF.spr[kind].stf2 * SQR(len);
		}
	}
  }
  
  // Extension of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		abpInd = abp.id[n];
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
		CONT(pArr[3] < 0);
		len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0), 0);
		lenEq = (abp.mId[n] == 0 && abp.mId[iAbp[pArr[3]]] == 0) 
				? motSA.cenDist : motSA.spr.eq;
		E[3] += 0.5 * motSA.spr.stf * SQR(len - lenEq);
	}
  }

  MPI_Gather(E, 4, MPI_DOUBLE, Eall, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    V4SET_ALL(Esum, 0.); 
    for(n = 0; n < nCpu; n++) {
        for(k = 0; k < 4; k++) {
            Esum[k] += Eall[4 * n + k];
        }
    }
	fOut = fopen(GenFileName("MechEall"), "a");
    fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayDouble(fOut, Esum, 4, 0);
    fclose(fOut);
    free(Eall);
  }
}

/*----------------------------- Energy calculation ---------------------------*/

/*------------------------ Internal stress measurement -----------------------*/

void RecordElasViscStressSubroutine1(double *rPnt, double *widDiv, int *dir, 
		double minDimDomH, int *indDomDiv, int *nSec, int mode) {
  int k;
  double dr[NDIM], ang, len;
  
  V3SET_ALL(indDomDiv, -1);
  if (mode == 0) {  
	FOR_NDIM(k) { 
		indDomDiv[k] = (int)((rPnt[k] - rGrid[k][0]) /  widDiv[k]);
		if (indDomDiv[k] < 0) { 
			indDomDiv[k] = (pbc[k] == 1) ? nSec[k] - 1 : 0;
		}
		else if (indDomDiv[k] >= nSec[k]) { 
			indDomDiv[k] = (pbc[k] == 1) ? 0 : nSec[k] - 1;
		}
	}
  }
  else {
	FOR_NDIM(k) {
		dr[k] = rPnt[k] - rGrid[k][0];
		CONT(k == dir[2]); 
		dr[k] -= dimDomH[k];
	}
	len = sqrt(SQR(dr[dir[0]]) + SQR(dr[dir[1]]));
	if (len >= minDimDomH) {
		return;
	}
	indDomDiv[0] = (int)(len / widDiv[0]);
	TrimIntVal(indDomDiv[0], 0, nSec[0] - 1);

	ang = atan(dr[dir[1]] / dr[dir[0]]);
	if (dr[dir[0]] < 0) { ang += PI; }
	if (ang < 0) { ang += 2 * PI; }
	indDomDiv[1] = (int)(ang / widDiv[1]);
	TrimIntVal(indDomDiv[1], 0, nSec[1] - 1);

	indDomDiv[2] = (int)(dr[dir[2]] / widDiv[2]);
	if (indDomDiv[2] < 0) { 
		indDomDiv[2] = (pbc[dir[2]] == 1) ? nSec[2] - 1 : 0;
	}
	else if (indDomDiv[2] >= nSec[2]) { 
		indDomDiv[2] = (pbc[dir[2]] == 1) ? 0 : nSec[2] - 1;
	} 
  }
}

void RecordElasViscStressSubroutine2(double stiff, double lenEq, 
		int *indDomDiv, double *rPnt, double *rPnt2, double *widDiv, int *dir, 
		double minDimDomH, int *nSec, int mode, int kind) {
  int n, k, dir2, CS, mm, l, idx, sft[NDIM];
  int ind[NDIM], indAdj[NDIM], indDomDivAdj[NDIM], begin[NDIM], end[NDIM];
  double len, len2, slope, ang, aa, bb, cc, dd, diff, low[NDIM], hi[NDIM];
  double dr[NDIM], dR[2][NDIM], rSec[NDIM], rAdj[2][NDIM], rCrs[NDIM];

  RecordElasViscStressSubroutine1(rPnt2, widDiv, dir, minDimDomH, 
		indDomDivAdj, nSec, mode);
  FOR_NDIM(k) {
	BREAK(indDomDiv[k] != indDomDivAdj[k]);
  }
  if (k == NDIM) { return; }
  FOR_NDIM(k) {
	rAdj[1][k] = rPnt2[k] - rGrid[k][0];
	rAdj[0][k] = rPnt[k] - rGrid[k][0];
	diff = rPnt2[k] - rPnt[k];
	if (pbc[k] == 1 && fabs(diff) > dimDomH[k]) {
		sft[k] = 1;
		if (diff >= dimDomH[k]) { rAdj[1][k] -= dimDom[k]; }
		else if (diff < -1. * dimDomH[k]) { rAdj[0][k] -= dimDom[k]; }
	}
	else { sft[k] = 0; }
	dr[k] = rAdj[1][k] - rAdj[0][k];
  }	

  FOR_NDIM(k) {
	if (((mode == 0 || (mode == 1 && k == dir[2])) && sft[k] == 1) 
			|| (mode == 1 && k == dir[1] && (double)abs(indDomDivAdj[k] 
			- indDomDiv[k]) >= (double)nSec[k] / 2.)) {
		begin[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDivAdj[k] 
				: indDomDiv[k];
		end[k] = ((indDomDiv[k] <= indDomDivAdj[k]) ? indDomDiv[k] 
				: indDomDivAdj[k]) + nSec[k];
	}
	else {
		begin[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDiv[k] 
				: indDomDivAdj[k];
		end[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDivAdj[k] 
				: indDomDiv[k];
	}
  }
  for(ind[0] = begin[0]; ind[0] <= end[0]; ind[0]++) {
	for(ind[1] = begin[1]; ind[1] <= end[1]; ind[1]++) {
		for(ind[2] = begin[2]; ind[2] <= end[2]; ind[2]++) {
			FOR_NDIM(k) {
				indAdj[k] = ind[k] % nSec[k];
			}
			V3IND_BACK_INT(idx, indAdj, nSec);
			V3MUL(rSec, indAdj, widDiv);
			if (mode == 0) {
				FOR_NDIM(k) {
					CONT(pbc[k] == 0 && indAdj[k] == 0);
					CONT((rAdj[0][k] - rSec[k]) * (rAdj[1][k] - rSec[k]) >= 0);
					CS = 1;
					for(n = 1; n < 3; n++) {
						dir2 = (k + n) % NDIM;
						rCrs[dir2] = (rSec[k] - rAdj[0][k]) * dr[dir2] / dr[k] 
								+ rAdj[0][dir2];
						if (rCrs[dir2] < rSec[dir2] || rCrs[dir2] >= rSec[dir2] 
								+ widDiv[dir2]) {
							CS = 0;
							break;
						}
					}
					CONT(CS == 0);
					rCrs[k] = rSec[k];
					len = V3LEN(dr);
					P2A(elasStre,idx,k,NDIM) += stiff * (len - lenEq) 
							* fabs(dr[k]) / len;
					P2A(elasStre2,idx, kind * NDIM + k,NDIM * 4) 
							+= stiff * (len - lenEq) * fabs(dr[k]) / len;

				}  
			}	
			else {
				FOR_NDIM(k) {
					for(n = 0; n < 2; n++) {
						dR[n][dir[k]] = rAdj[n][dir[k]] 
								- (k == dir[2] ? 0 : dimDomH[dir[k]]);
					}
					if (dR[0][dir[k]] < dR[1][dir[k]]) {
						low[k] = dR[0][dir[k]];
						hi[k] = dR[1][dir[k]];
					}
					else {
						low[k] = dR[1][dir[k]];
						hi[k] = dR[0][dir[k]];
					}
				}
				dd = dR[1][dir[0]] * dR[0][dir[1]] 
						- dR[0][dir[0]] * dR[1][dir[1]];
				FOR_NDIM(k) {
					CS = 1;
					if (k == 0) {
						len = sqrt(SQR(dR[0][dir[0]]) + SQR(dR[0][dir[1]]));
						len2 = sqrt(SQR(dR[1][dir[0]]) + SQR(dR[1][dir[1]]));
						CONT((len - rSec[k]) * (len2 - rSec[k]) >= 0);
						aa = 1 + SQR(dr[dir[1]] / dr[dir[0]]);
						bb = dr[dir[1]] * dd / SQR(dr[dir[0]]);
						cc = SQR(dd / dr[dir[0]]) - SQR(rSec[k]);
						for(n = 0; n < 2; n++) {
							rCrs[dir[0]] = (-1. * bb + (n == 0 ? 1. : -1.) 
									* sqrt(bb * bb - aa * cc)) / aa;
							BREAK(rCrs[dir[0]] >= low[dir[0]] 
									&& rCrs[dir[0]] < hi[dir[0]]);
						}
						for(n = 1; n < 3; n++) {
							rCrs[dir[n]] = dr[dir[n]] / dr[dir[0]] 
									* (rCrs[dir[0]] - dR[0][dir[0]]) 
									+ dR[0][dir[n]];
						}
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						if (ang < rSec[1] || ang >= rSec[1] + widDiv[1]
								|| rCrs[dir[2]] < rSec[2] 
								|| rCrs[dir[2]] >= rSec[2] + widDiv[2]) {
							CS = 0;
						}
					}
					else if (k == 1) {
						if (fabs(rSec[1] - 0.5 * PI) < 0.01 
								|| fabs(rSec[1] - 1.5 * PI) < 0.01) {
							CONT(fabs(dr[dir[0]]) < 0.01);
							CONT(dR[0][dir[0]] * dR[1][dir[0]] >= 0.);
							rCrs[dir[0]] = 0.;
						}
						else {
							slope = tan(rSec[1]);
							CONT(fabs(slope - dr[dir[1]] / dr[dir[0]]) 
									< 0.01);
							rCrs[dir[0]] = dd / (slope * dr[dir[0]] 
									- dr[dir[1]]);
							CONT(rCrs[dir[0]] < low[dir[0]] 
									|| rCrs[dir[0]] >= hi[dir[0]]);
						}
						for(n = 1; n < 3; n++) {
							rCrs[dir[n]] = dr[dir[n]] / dr[dir[0]] 
									* (rCrs[dir[0]] - dR[0][dir[0]]) 
									+ dR[0][dir[n]];
						}
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						CONT(ang < rSec[1] || ang >= rSec[1] + widDiv[1]);

						len = sqrt(SQR(rCrs[dir[0]]) + SQR(rCrs[dir[1]]));
						if (len < rSec[0] || len >= rSec[0] + widDiv[0]
								|| rCrs[dir[2]] < rSec[2] 
								|| rCrs[dir[2]] >= rSec[2] + widDiv[2]) {
							CS = 0;
						}
					}
					else {
						CONT(pbc[k] == 0 && indAdj[k] == 0);
						CONT((rAdj[0][k] - rSec[k]) 
								* (rAdj[1][k] - rSec[k]) >= 0);
						CS = 1;
						rCrs[dir[0]] = (rSec[k] - rAdj[0][k]) * dr[dir[0]] 
								/ dr[k] + dR[0][dir[0]];
						rCrs[dir[1]] = (rSec[k] - rAdj[0][k]) * dr[dir[1]] 
								/ dr[k] + dR[0][dir[1]];
						len = sqrt(SQR(rCrs[dir[0]]) + SQR(rCrs[dir[1]]));
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						if (len < rSec[0] || len >= rSec[0] + widDiv[0] 
								|| ang < rSec[1] || ang >= rSec[1] 
								+ widDiv[1]) {
							CS = 0;
						}
					}
					CONT(CS == 0);
					len = V3LEN(dr);
					P2A(elasStre,idx,k,NDIM) += stiff * (len - lenEq) 
							* fabs(dr[k]) / len;
					P2A(elasStre2,idx, kind * NDIM + k,NDIM * 4) 
							+= stiff * (len - lenEq) * fabs(dr[k]) / len;
				}  
			}
		}
	}
  }
}

// mode = 0: Cartesian, 1: cylindrical
void RecordElasViscStress(int period, int *nSec, int mode) {
  int mm, m, n, k, end, nSecP, nSecPN, CS, side;
  int abpInd, *pArr, ind, indDomDiv[NDIM], dir[NDIM];
  double fac, minDimDomH, widDiv[NDIM];
  double dr[NDIM], dr2[NDIM], drAxis[NDIM], rPos[NDIM], rPos2[NDIM];
  double *rPnt, *fPnt, *rPntAdj, *rPntAdjAdj, *streAll, *streSum, *sendArr;
  FILE *fOut;

  minDimDomH = POS_LARGE_VALUE;
  nSecP = V3PROD(nSec);
  if (mode == 0) {
	V3DIV(widDiv, dimDom, nSec);
  }
  else {
	V3SET(dir, 0, 1, 2);
	FOR_NDIM(k) {
		CONT(!(k != dir[2] && dimDomH[k] < minDimDomH));
		minDimDomH = dimDomH[k];
	}
	widDiv[0] = minDimDomH / nSec[0];
	widDiv[1] = 2 * PI / nSec[1];
	widDiv[2] = dimDom[dir[2]] / nSec[2];
  }

  for(n = 0; n < nActMe + nAbpMe; n++) {
	CS = 1;
	if (n < nActMe) {
		CONT(ISACTM(n));
		fPnt = &P2(act.f,n,0);
		rPnt = &P2(act.r,n,0);
		if (act.fix[n] > -1 && bndReb.gTglSpr == 0) { CS = -1; }
	}
	else {
		abpInd = n - nActMe;
		CONT(ISABPM(abpInd));
		pArr = &P2A(abp.ch,abpInd,0,nChAb);
		fPnt = &P2(abp.f,abpInd,0);
		rPnt = &P2(abp.r,abpInd,0);
	}
	RecordElasViscStressSubroutine1(rPnt, widDiv, dir, minDimDomH, 
			indDomDiv, nSec, mode);
	V3IND_BACK_INT(ind, indDomDiv, nSec);
	CONT(indDomDiv[0] < 0);
	// Viscous
	if (CS == 1 && isnan(fPnt[0]) == 0) {
		VV3ADD(&P2A(viscStre,ind,0,NDIM), fPnt);
	}
	// Elastic
	if (n < nActMe) {
		end = 1;
	}
	else {
		end = (motSA.gTgl != 0) ? nChAb - 1 : 2;
	}
	for(m = 0; m < end; m++) {
		if (n < nActMe) {
			CONT(P2A(act.ch,n,0,nChAc) < 0);
			rPntAdj = &P2(act.r,iAct[P2A(act.ch,n,0,nChAc)],0);
		}
		else {
			CONT(m == 2);
			CONT(pArr[m] < 0);
			CONT(pArr[2] != 2 && m > 2);
			rPntAdj = (m < 2) ? &P2(act.r,iAct[pArr[m]],0) 
					: &P2(abp.r,iAbp[pArr[m]],0);
		}
		// Motor arm
		if (n >= nActMe && m < 2) {
			side = FindAbpActinChain(iAct[pArr[m]], abp.id[abpInd], 0);
			rPntAdjAdj = &P2(act.r,
					iAct[P2A(act.ch,iAct[pArr[m]],0,nChAc)],0);
			CalcPosOnActSegSide(rPntAdj, rPntAdjAdj, rPos, side);
  			CalcVec(dr, rPos, rPnt);

			if (ISMTF(pArr[2])) {
				CalcVec(drAxis, rPntAdj, rPntAdjAdj);
				fac = V3DOT(dr, drAxis) / V3LEN_SQ(drAxis);
				VS3SUB(rPos2, rPos, drAxis, fac);
				ApplyBoundCondVector(rPos2, -1, 0);
			}
			RecordElasViscStressSubroutine2(abpF.spr[pArr[2]].stf,  
					abpF.spr[pArr[2]].eq, indDomDiv, rPnt, ((ISMTF(pArr[2]))
					? rPos2 : rPos), widDiv, dir, minDimDomH, nSec, mode,
					(pArr[2] == 2) ? 2 : 1);
			if (ISMTF(pArr[2])) {
				RecordElasViscStressSubroutine2(abpF.spr[pArr[2]].stf2,  
						0., indDomDiv, rPnt, rPos, widDiv, dir, minDimDomH,
						nSec, mode, 2);
			}
		}
		else {
			// F-actin
			if (n < nActMe) { 
				RecordElasViscStressSubroutine2(actF.spr.stf, actF.spr.eq, 
						indDomDiv, rPnt, rPntAdj, widDiv, dir, minDimDomH, 
						nSec, mode, 0);
			}
			// Motor backbone
			else {
				RecordElasViscStressSubroutine2(motSA.spr.stf, motSA.spr.eq,
						indDomDiv, rPnt, rPntAdj, widDiv, dir, minDimDomH, 
						nSec, mode, 3);
			}
		}
	}
  }
  nSecPN = nSecP * NDIM;
  if (currTimeStep % period == 0) {
    MALLOC(sendArr,double,nSecPN * (2 + 4));
    if (rank == 0) {
        MALLOC(streAll,double,nSecPN * (2 + 4) * nCpu);
        MALLOC(streSum,double,nSecPN * (2 + 4));
    }
    for(n = 0; n < nSecPN; n++) {
        sendArr[n] = elasStre[n];
        sendArr[n + nSecPN] = viscStre[n];
    }
    for(n = 0; n < nSecPN * 4; n++) {
		sendArr[n + nSecPN * 2] = elasStre2[n];
	}
	
    MPI_Gather(sendArr, nSecPN * (2 + 4), MPI_DOUBLE, streAll,
            nSecPN * (2 + 4), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for(n = 0; n < nSecPN * (2 + 4); n++) {
            streSum[n] = 0.;
            for(k = 0; k < nCpu; k++) {
                streSum[n] += streAll[k * nSecPN * (2 + 4) + n];
            }
            streSum[n] /= (double)period;
            streSum[n] = F_S2PN(streSum[n]);
        }
        fOut = fopen(GenFileName("ElasViscStre"), "a");
        if (currTimeStep / period == 1) {
            fprintf(fOut, "%d\t", nSecP * 2);
            Fprintf1dArrayIntWoRet(fOut, nSec, NDIM);
            FOR_NDIM(k) {
                fprintf(fOut, "%g\t", ((mode == 1 && k == 1)
                        ? RAD2DEG(widDiv[k]) : L_S2UM(widDiv[k])));
            }
            fprintf(fOut, "\n");
        }
        for(n = 0; n < nSecP * 2; n++) {
            V3IND_ASSIGN_INT((n % nSecP), nSec, 1, indDomDiv);
            fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\n", n, indDomDiv[0],
                    indDomDiv[1], indDomDiv[2], P2(streSum,n,0),
                    P2(streSum,n,1), P2(streSum,n,2));
        }
        fclose(fOut);

        fOut = fopen(GenFileName("ElasForce"), "a");
        if (currTimeStep / period == 1) {
            fprintf(fOut, "%d\t", nSecP);
            Fprintf1dArrayIntWoRet(fOut, nSec, NDIM);
            FOR_NDIM(k) {
                fprintf(fOut, "%g\t", ((mode == 1 && k == 1)
                        ? RAD2DEG(widDiv[k]) : L_S2UM(widDiv[k])));
            }
			for(k = 0; k < NDIM * 4 - 3; k++) {
				fprintf(fOut, "0\t");
			}
            fprintf(fOut, "\n");
        }
        for(n = 0; n < nSecP; n++) {
            V3IND_ASSIGN_INT(n, nSec, 1, indDomDiv);
            fprintf(fOut, "%d\t%d\t%d\t%d\t", n, indDomDiv[0],
                    indDomDiv[1], indDomDiv[2]);
			for(k = 0; k < NDIM * 4; k++) {
            	fprintf(fOut, "%g\t", streSum[n * NDIM * 4 + k + nSecPN * 2]);
			}
			fprintf(fOut, "\n");
        }
        fclose(fOut);

        free(streSum);
        free(streAll);
    }
    free(sendArr);
	SetAllValue1dArrayDouble(viscStre, nSecPN, 0.);
	SetAllValue1dArrayDouble(elasStre, nSecPN, 0.);
	SetAllValue1dArrayDouble(elasStre2, nSecPN * 4, 0.);
  }
}

/*--------------------------- Accumulated information ------------------------*/

// Accumulate forces acting on actin and ABP. This information is not recorded 
// to file, but used for coloring methods via VMD.
void RecordAccuLengthForces(void) {
  int n, k, cntCh;
  int actInd, nextActInd, nextAbpInd;
  double lenAccu, len;
  // actin
  FOR_ACTME(n) {
	CONT(ISACTM(n));
	VV3ADD(&P2A(recAct.allF,n,0,NDIM + 1), &P2(act.f,n,0));
	P2A(recAct.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(act.f,n,0));
	nextActInd = P2A(act.ch,n,0,nChAc);
	if (nextActInd > -1) {
		len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[nextActInd],0), 0);
		recAct.len[n] += len;
		recAct.cnt[n] += 1;
	}
	else {
		recAct.len[n] = 0.;
		recAct.cnt[n] = 0;
	}
  } 
  // ABP
  FOR_ABPME(n) {
	CONT(ISABPM(n));
	VV3ADD(&P2A(recAbp.allF,n,0,NDIM + 1), &P2(abp.f,n,0));
	P2A(recAbp.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(abp.f,n,0));
	cntCh = 0;
	lenAccu = 0.;
	for(k = 0; k < 2; k++) {
		actInd = P2A(abp.ch,n,k,nChAb);
		CONT(actInd < 0);
		len = CalcDistActinAbp(actInd, abp.id[n], 0);
		lenAccu += len;
		cntCh++;
	}
	if (cntCh > 0) {
		P2A(recAbp.len,n,0,recAbp.nL) += lenAccu / (double)cntCh;
		recAbp.cnt[n] += 1;
	}
	else {
		P2A(recAbp.len,n,0,recAbp.nL) = 0.;
		if (!(ISMTF(K_ABP(n)))) {
			recAbp.cnt[n] = 0;
		}
	}
	CONT(!(ISMTF(K_ABP(n))));
	nextAbpInd = P2A(abp.ch,n,3,nChAb);
	if (nextAbpInd > -1) {
		len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[nextAbpInd],0), 0);
		P2A(recAbp.len,n,1,recAbp.nL) += len;
		if (cntCh == 0) { 
			recAbp.cnt[n] += 1;
		}
	}
	else {
		P2A(recAbp.len,n,1,recAbp.nL) = 0.;
		if (cntCh == 0) {
			recAbp.cnt[n] = 0;
		}
	}
  }
}

void ResetAccuLengthForces(void) {
  int n;

  FOR_ACTC(n) {
	recAct.sprF[n] = 0.;
	recAct.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAct.allF,n,0,NDIM + 1), 0.);
	recAct.len[n] = 0.;
	recAct.cnt[n] = 0.;
  }
  FOR_ABPC(n) {
	recAbp.sprF[n] = 0.;
	recAbp.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAbp.allF,n,0,NDIM + 1), 0.);
	SetAllValue1dArrayDouble(&P2A(recAbp.len,n,0,recAbp.nL), recAbp.nL, 0.);
	recAbp.cnt[n] = 0.;
  }
}

// Calculate and record the longitudinal forces acting on ACP or motor. This
// information is very useful to analyze the activity of motors since 
// the longitudinal forces are known to affect the rate of walking events.
void RecordAccuLongSpringForces(int period) {
  int n, k, *nAbpAll, *pArr, locActInd, side, nR;
  double *f, *fAbp, *rPnt[2], *rBnd, *rBndAll;
  double rPos[NDIM], rPos2[NDIM], dr[2][NDIM], fac;
  FILE *fOut;

  nR = (motSA.gTgl != 0 && nMot > 0) ? 4 : 2;

  // Accmulate the calculated longitudinal forces.
  FOR_ABPME(n) {
	CONT(ISABPM(n));
    for (k = 0; k < 2; k++) {
		CONT(P2A(abp.ch,n,k,nChAb) < 0);
        // find force
        f = &P2A(recInstSprFabp,n,k * 2,4);
		P2A(recLongSprFabp,n * 2 + k,0,2) += f[0];
		P2A(recLongSprFabp,n * 2 + k,1,2) += f[1];
	}
  }
  if (currTimeStep % period == 0) {
	MALLOC(nAbpAll,int,nCpu);
	if (rank == 0) { 
		MALLOC(fAbp,double,nAbp * 2 * 2); 
		MALLOC(chAbp,int,nAbp * nChAb); 
		MALLOC(rAbp,double,nAbp * NDIM); 
		MALLOC(rBndAll,double,nAbp * NDIM * nR);
		SetAllValue1dArrayDouble(rBndAll, nAbp * NDIM * nR, 0.);
	}
	MALLOC(rBnd, double,nAbpMe * NDIM * nR);
	SetAllValue1dArrayDouble(rBnd, nAbpMe * NDIM * nR, 0.);
	
	FOR_ABPME(n) {
		CONT(ISABPM(n));
		pArr = &P2A(abp.ch,n,0,nChAb);
    	for (k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			locActInd = iAct[pArr[k]];
			side = FindAbpActinChain(locActInd, abp.id[n], 0);
			rPnt[0] = &P2(act.r,locActInd,0);
			rPnt[1] = &P2(act.r,iAct[P2A(act.ch,locActInd,0,nChAc)],0);
			CalcPosOnActSegSide(rPnt[0], rPnt[1], rPos, side);
	  		CalcVec(dr[0], rPos, &P2(abp.r,n,0));
			if (ISMTF(pArr[2])) {
				CalcVec(dr[1], rPnt[0], rPnt[1]);
				fac = V3DOT(dr[0], dr[1]) / V3LEN_SQ(dr[1]);
				VS3SUB(rPos2, rPos, dr[1], fac);
				ApplyBoundCondVector(rPos2, -1, 0);
			}
			V3COPY(&P2A(rBnd,n,(k * nR / 2) * NDIM,NDIM * nR), rPos);	
			if (ISMTF(pArr[2])) {
				V3COPY(&P2A(rBnd,n,(k * nR / 2 + 1) * NDIM,NDIM * nR), rPos2);
			}
		}
  	}

	MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Gather2dArrayDouble(nAbpMe, nAbpAll, 4, abp.id, recLongSprFabp, fAbp);	
	Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);
	Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM, abp.id, abp.r, rAbp);
	Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM * nR, abp.id, rBnd, rBndAll);
    for(n = 0; n < nAbpC * 4; n++) {
        recLongSprFabp[n] = 0.;
    }

	if (rank == 0) {
		for(n = 0; n < nAbp * 4; n++) {
			fAbp[n] = F_S2PN(fAbp[n]) / (double)period;
		}
		fOut = fopen(GenFileName("AbpLongFor"), "a");
		if (currTimeStep / period == 1) {
    		Fprintf1dArrayIntWFil(fOut, &nAbp, 1, 0, 19 + nR);
		}
		
		FOR_ABP(n) {
			fprintf(fOut, "%d\t", n);
			for(k = 0; k < 4; k++) {
				fprintf(fOut, "%g\t", P2A(fAbp,n,k,4));
			}
			for(k = 0; k < 3; k++) {
				fprintf(fOut, "%d\t", P2A(chAbp,n,k,nChAb));
			}
			FOR_NDIM(k) {
				fprintf(fOut, "%g\t", P2(rAbp,n,k));
			}	
			for(k = 0; k < NDIM * nR; k++) {
				fprintf(fOut, "%g\t", P2A(rBndAll,n,k,NDIM * nR));
			}	
			fprintf(fOut, "\n");
		}

		fclose(fOut);
		free(nAbpAll);
		free(fAbp);
		free(chAbp);
		free(rAbp);
		free(rBndAll);
	}
	free(rBnd);
  }
}

/*--------------------------- Accumulated information ------------------------*/

/*------------------ Information for filaments and segments  -----------------*/

void RecordIndvSegFilaInformation(int period) {
  int n, k, l, curr, ind, iPercSum, iSuppSum, cnt, cntFor, CS, side;
  int abpInd, actInd, prevAbpInd, prevActInd, sft[NDIM], cntXlink[3];
  int *iFila, *nActAll, *nAbpAll, *cntAct, *cntAbp, *iMot;
  int *binAct, *binAbp, *binActAll, *binAbpAll, *pArr, *pArr2;
  double len, lenCtr, lenEE, lenSeg, ang, angSum[NDIM], crsAng[2];
  double tensFactSum, tensFabpSum, bendFactSum, bendFabpSum;
  double dr[NDIM], dr2[NDIM], drCrs[NDIM], drSum[NDIM];
  double viscFsum[NDIM], rPos[2][NDIM], ratio[2];
  double *fAct, *fAbp, *fBendAct, *fBendAbp;
  double *lAct, *lAbp, *tensF, *bendF, *viscF;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(binAct,int,nActC);
  MALLOC(binAbp,int,nAbpC);
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	MALLOC(fAct,double,nAct * (NDIM + 1));
	MALLOC(fAbp,double,nAbp * (NDIM + 1));
	MALLOC(fBendAct,double,nAct);
	MALLOC(fBendAbp,double,nAbp);
	MALLOC(rAct,double,nAct * NDIM);
	MALLOC(rAbp,double,nAbp * NDIM);
	MALLOC(lAct,double,nAct);
	MALLOC(lAbp,double,nAbp * recAbp.nL);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(cntAct,int,nAct);
	MALLOC(cntAbp,int,nAbp);
	MALLOC(tensF,double,nAct + nAbp);
	MALLOC(bendF,double,nAct + nAbp);
	MALLOC(viscF,double,(nAct + nAbp) * NDIM);
	MALLOC(fixAct,int,nAct);
	MALLOC(iFila,int,nAct);
	MALLOC(binActAll,int,nAct);
	MALLOC(binAbpAll,int,nAbp);
	if (motSA.gTgl != 0) {
		MALLOC(iMot,int,nAbp);
	}
  }
  FOR_ACTME(n) { 
	binAct[n] = BinaryPackActinChainArray(act.id[n]);
  }
  FOR_ABPME(n) { 
	binAbp[n] = BinaryPackAbpChainArray(abp.id[n]);
  }
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayDouble(nActMe, nActAll, NDIM + 1, act.id, recAct.allF, fAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM + 1, abp.id, recAbp.allF, fAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.bendF, fBendAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, 1, abp.id, recAbp.bendF, fBendAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.len, lAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, recAbp.nL, abp.id, recAbp.len, lAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, recAct.cnt, cntAct);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, recAbp.cnt, cntAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, binAct, binActAll);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, binAbp, binAbpAll);

  if (rank == 0) {
	if (motSA.gTgl != 0) {
		memset(iMot, -1, sizeof(int) * nAbp);
		cnt = 0;
		FOR_ABP(n) {
			pArr = &P2A(chAbp,n,0,nChAb);
			CONT(!(pArr[3] > -1 && pArr[4] < 0));
			curr = n;
			while(P2A(chAbp,curr,3,nChAb) > -1) {
				iMot[curr] = cnt;
				curr = P2A(chAbp,curr,3,nChAb);
			}
			iMot[curr] = cnt;
			cnt++;
		}
	}
	fOut = fopen(GenFileName("InfoIndv"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t%d\t%d\t", nAct + nAbp, nAct, nAbp);
		Fprintf1dFillerInt(fOut, 0, 16, 0);
	}
	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) {
			fprintf(fOut, "%d\t", n);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 7, 0);
			tensF[n] =	0.;
			bendF[n] =	0.;
			V3SET_ALL(&P2(viscF,n,0), 0.);
			continue;
		}
		else if (cntAct[n] == 0) {
			tensF[n] =	0.;
			bendF[n] =	0.;
			V3SET_ALL(&P2(viscF,n,0), 0.);
		}
		else {
			tensF[n] = F_S2PN(actF.spr.stf * (lAct[n] / (double)cntAct[n] 
					- actF.spr.eq));
			bendF[n] = F_S2PN(fBendAct[n] / (double)cntAct[n]);
			VS3COPY(&P2(viscF,n,0), &P2A(fAct,n,0,NDIM + 1), INV(cntAct[n]));
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			CalcUnitVec(dr, &P2(rAct,n,0), &P2(rAct,pArr[k],0));
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);

		FOR_NDIM(k) { P2(viscF,n,k) = F_S2PN(P2(viscF,n,k)); }
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n, P2(rAct,n,0), 
				P2(rAct,n,1), P2(rAct,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[n], bendF[n], P2(viscF,n,0), P2(viscF,n,1), P2(viscF,n,2),
				pArr[0], pArr[1], iPerc[n], iSupp[n], fixAct[n], iFila[n], 
				binActAll[n]);
	}
	FOR_ABP(n) {
		ind = n + nAct;
		pArr = &P2A(chAbp,n,0,nChAb);
		if ((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0)) {
			fprintf(fOut, "%d\t", ind);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 7, 0);
			tensF[ind] = 0.;
			bendF[ind] = 0.;
			V3SET_ALL(&P2(viscF,ind,0), 0.);
			continue;
		}
 		else if (cntAbp[n] == 0) {
			tensF[ind] = 0.;
			bendF[ind] = 0.;
			V3SET_ALL(&P2(viscF,ind,0), 0.);
		}
		else {
			if (ISMTF(pArr[2])) {
				tensF[ind] = F_S2PN(motSA.spr.stf 
						* (P2A(lAbp,n,1,recAbp.nL) / (double)cntAbp[n] 
						- motSA.spr.eq));
			}
			else {
				tensF[ind] = F_S2PN(abpF.spr[pArr[2]].stf 
						* (P2A(lAbp,n,0,recAbp.nL) / (double)cntAbp[n]
						- abpF.spr[pArr[2]].eq));
			}
			bendF[ind] = F_S2PN(fBendAbp[n] / (double)cntAbp[n]);
			VS3COPY(&P2(viscF,ind,0), &P2A(fAbp,n,0,NDIM + 1), INV(cntAbp[n]));
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		ind = (ISMTF(pArr[2])) ? 3 : 0;
		for(k = 0; k < 2; k++) {
			CONT(pArr[k + ind] < 0);
			if (ind == 0) {
				CalcUnitVecActinAbp(dr, pArr[k], n, 1);
			}
			else {
				CalcUnitVec(dr, &P2(rAbp,pArr[k + ind],0), &P2(rAbp,n,0)); 
			}
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);
		ind = n + nAct;
		FOR_NDIM(k) { P2(viscF,ind,k) = F_S2PN(P2(viscF,ind,k)); }
		for(k = 0; k < 2; k++) {
			side = FindAbpActinChain(pArr[k], n, 1);
			ratio[k] = (side > -1) ? (double)((side - 2) / nChAcY) 
					/ (double)nChAcX : -1.;
		}
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%g\t%g\t%d\n", ind, P2(rAbp,n,0), 
				P2(rAbp,n,1), P2(rAbp,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[ind], bendF[ind], P2(viscF,ind,0), P2(viscF,ind,1), 
				P2(viscF,ind,2), pArr[0], pArr[1], iPerc[ind], iSupp[ind], 
				ratio[0], ratio[1],	binAbpAll[n]);
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFila"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		cnt++;
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 26);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		tensFactSum = 0.;
		tensFabpSum = 0.;
		bendFactSum = 0.;
		bendFabpSum = 0.;
		lenSeg = 0.;
		lenCtr = 0.;
		iPercSum = 0;
		iSuppSum = 0;
		V3SET_ALL(angSum, 0.);
		V3SET_ALL(cntXlink, 0);
		V3SET_ALL(drSum, 0.);
		V3SET_ALL(viscFsum, 0.);
		V3SET_ALL(sft, 0);
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			lenCtr += len;
			CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
			if (pArr[1] > -1) {
				// Calculate bending angles
				CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
				V3CROSS(drCrs, dr, dr2);
				NormVec(drCrs);
				ang = V3ANG(dr, dr2);
				VVS3ADD(angSum, drCrs, ang);
			}
			VVS3ADD(drSum, dr, INV(len));
			tensFactSum += tensF[curr];
			bendFactSum += bendF[curr];
			for(k = 2; k < nChAc; k++) {
				abpInd = pArr[k];
				CONT(abpInd < 0);
				pArr2 = &P2A(chAbp,abpInd,0,nChAb);
				CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
				tensFabpSum += tensF[nAct + abpInd];
				bendFabpSum += bendF[nAct + abpInd];
				iSuppSum += (iSupp[nAct + abpInd] > -1) ? 1 : 0;
				cntXlink[(int)(pArr2[2] / 2)]++; 
				CONT(ISMTF(pArr2[2]));
				actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
				CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
						&P2(rAct,actInd,0));
				if (V3DOT(dr, dr2) > 0.) { cntXlink[2]++; }
			}
			VV3ADD(viscFsum, &P2(viscF,curr,0));
			iPercSum += iPerc[curr] % 10;
			iSuppSum += (iSupp[curr] > -1) ? 1 : 0;
			lenSeg += 1.;
			curr = pArr[0];
		}
		// Normalize the orientation vector.
		NormVec(drSum);
		// Measure the end-to-end distance of the filament.
		V3SUB(dr, &P2(rAct,curr,0), &P2(rAct,n,0));
		VSV3ADD(dr, dr, sft, dimDom);
		lenEE = V3LEN(dr);
		if (lenSeg > 0.) {
			tensFactSum /= lenSeg;
			bendFactSum /= lenSeg;
			V3SCALE(viscFsum, INV(lenSeg));
		}
		if (cntXlink[0] + cntXlink[1] > 0) {
			tensFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
			bendFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
		}
		fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
				"\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\t%d\n", cnt, 
				iFila[n], n, curr, drSum[0], drSum[1], drSum[2], tensFactSum, 
				tensFabpSum, bendFactSum, bendFabpSum, L_S2UM(lenSeg), 
				L_S2UM(lenCtr), L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
				viscFsum[0], viscFsum[1], viscFsum[2], cntXlink[0], cntXlink[1],
				cntXlink[2], iPercSum, iPerc[n] / 10, iSuppSum);
		cnt++;
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFilaSeg"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = P2A(chAct,curr,2 + k * nChAcY + l,nChAc);
					CONT(abpInd < 0);
					CONT(!(ISMTF(P2A(chAbp,abpInd,2,nChAb)))
							&& P2A(chAbp,abpInd,1,nChAb) < 0);
					if (CS == 1) { 
						cnt++;
					}
					CS = 1;
					prevAbpInd = abpInd;
				}
			}
			curr = P2A(chAct,curr,0,nChAc);
		}
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 22);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			if (CS == 1) {
				VVS3ADD(drSum, dr, INV(len));
				tensFactSum += tensF[curr];
				bendFactSum += bendF[curr];
				VV3ADD(viscFsum, &P2(viscF,curr,0));
				cntFor++;
				iPercSum += iPerc[curr] % 10;
				iSuppSum += (iSupp[curr] > -1) ? 1 : 0;
			}
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = pArr[2 + k * nChAcY + l];
					CONT(abpInd < 0);
					pArr2 = &P2A(chAbp,abpInd,0,nChAb);
					CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
					CalcPosOnActSeg(&P2(rAct,curr,0), &P2(rAct,pArr[0],0), 
							rPos[1], (double)k / (double)nChAcX);
					crsAng[1] = -1.;
					if (!(ISMTF(pArr2[2]))) {
						actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
						CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
								&P2(rAct,actInd,0));
						crsAng[1] = V3ANG(dr, dr2);
					}
					if (CS == 1) { 
						if (curr != prevActInd) {
							CheckCrossBound(sft, &P2(rAct,curr,0), rPos[1]);
						}
						V3SUB(dr2, rPos[1], rPos[0]);
						VSV3ADD(dr2, dr2, sft, dimDom);
						lenEE = V3LEN(dr2);
						if (cntFor > 0) {
							tensFactSum /= (double)cntFor;
							bendFactSum /= (double)cntFor;
							V3SCALE(viscFsum, INV((double)cntFor));
						}
						NormVec(drSum);
						fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t"
								"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\t"
								"%g\n", cnt, iFila[n], prevActInd, curr, 
								drSum[0], drSum[1], drSum[2], tensFactSum, 
								bendFactSum, L_S2UM(lenSeg), L_S2UM(lenCtr), 
								L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
								viscFsum[0], viscFsum[1], viscFsum[2],
								prevAbpInd * 10 + P2A(chAbp,prevAbpInd,2,nChAb),
								abpInd * 10 + pArr2[2], crsAng[0], crsAng[1]);
						cnt++;
					}
					lenCtr = 0.;
					lenSeg = 0.;
					tensFactSum = tensF[curr];
					bendFactSum = bendF[curr];
					V3COPY(viscFsum, &P2(viscF,curr,0));
					cntFor = 1;
					V3SET_ALL(angSum, 0.);
					VS3COPY(drSum, dr, INV(len));
					iPercSum = iPerc[curr] % 10;
					iSuppSum = (iSupp[curr] > -1) ? 1 : 0;
					V3SET_ALL(sft, 0);
					prevAbpInd = abpInd;
					prevActInd = curr;
					V3COPY(rPos[0], rPos[1]);
					crsAng[0] = crsAng[1];
					CS = 1;
				}
				if (CS == 1) {
					lenCtr += len / (double)nChAcX;
					lenSeg += 1. / (double)nChAcX;
				}
			}
			if (CS == 1) { 
				CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
				if (pArr[1] > -1) {
					// Calculate bending angles
					CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
					V3CROSS(drCrs, dr, dr2);
					NormVec(drCrs);
					ang = V3ANG(dr, dr2);
					VVS3ADD(angSum, drCrs, ang);
				}
			}
			curr = pArr[0];
		}
	}
	fclose(fOut);
  }

  if (rank == 0) {
	free(rAct);
	free(rAbp);
	free(fAct);
	free(fAbp);
	free(fBendAct);
	free(fBendAbp);
	free(lAct);
	free(lAbp);
	free(chAct);
	free(chAbp);
	free(cntAct);
	free(cntAbp);
	free(tensF);
	free(bendF);
	free(viscF);
	free(fixAct);
	free(iFila);
	free(binActAll);
	free(binAbpAll);
	if (motSA.gTgl != 0) {
		free(iMot);
	}
  }
  free(nActAll);
  free(nAbpAll);
  free(binAct);
  free(binAbp);
}

/*------------------ Information for filaments and segments  -----------------*/


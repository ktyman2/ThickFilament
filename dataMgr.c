// This file contains functions performing the modification or loading of 
// data.

/*------------------------ Related to network data ---------------------------*/

void ExtractConfigSubroutine(double *r, double *cBnd, double cNeiEdge, 
		int m, int *CS, int *CS2) {
  int k, CS3;

  ConvertRectDomainVector(r, 1);
  FOR_NDIM(k) {
	CS3 = 1;
	if (pbc[k] == 1 || (pbc[k] != 1 && iCell[k] > 0 
			&& iCell[k] < nCell[k] - 1)) {
		if (r[k] >= P2A(bnd.r,0,k,NDIM) && r[k] < P2A(bnd.r,1,k,NDIM)) { 
			(*CS)++; 
			CS3 = 0;
		}
	}
	else if ((iCell[k] == 0 && r[k] < P2A(bnd.r,1,k,NDIM)) 
			|| (iCell[k] == nCell[k] - 1 && r[k] >= P2A(bnd.r,0,k,NDIM))) {
		(*CS)++;
		CS3 = 0;
	}
	if (m == 1 && CS3 == 1 && nCell[k] > 1 
			&& ((r[k] >= P2A(cBnd,0,k,NDIM) - cNeiEdge
			&& r[k] < P2A(cBnd,0,k,NDIM)) || (r[k] >= P2A(cBnd,1,k,NDIM)  
			&& r[k] < P2A(cBnd,1,k,NDIM) + cNeiEdge))) {
		(*CS2)++;  
	}
  }
}

// After the data of a whole network are loaded at LoadConfig(), information
// of particles belonging to the current subdomain is extracted from the data
// here. The rest of data are discarded after this function.
void ExtractConfig(void) {
  int m, n, k, CS, CS2, ind, curr, nIFilaP, nActM, nAbpM, cnt, *pArr;
  double rBnd[NDIM*2], cNeiEdge, len;
  ListInt iFila;


  MALLOC(iFila.l, int, nAct);
  memset(iAct, -1, sizeof(int) * nAct);
  memset(iAbp, -1, sizeof(int) * nAbp);
  cNeiEdge = (motWalk.gTgl != 0) ? neiEdge + 1. : neiEdge;
  V6COPY(rBnd,bnd.r);
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 1));
	if (iCell[k] == 0 && nCell[k] > 1) 
	{ P2A(rBnd,0,k,NDIM) = rGrid[k][nGrid[k] - 1]; }
	else if (iCell[k] == nCell[k] - 1 && nCell[k] > 1) 
	{ P2A(rBnd,1,k,NDIM) = rGrid[k][0]; }
  }

  nActM = 0;
  memset(iFila.l, -1, sizeof(int) * nAct);
  iFila.c = 0;
  FOR_ACT(n) {
	pArr = &P2A(chAct,n,0,nChAc);
    // Count the number of free actins
    if (pArr[0] < 0 && pArr[1] < 0) {
        nActM++;
    }
    // "iFila.l" is the index of actin filaments.
    if (pArr[1] < 0 && pArr[0] > -1) {
        curr = n;
        while(P2A(chAct,curr,0,nChAc) > -1) {
            iFila.l[curr] = iFila.c;
            curr = P2A(chAct,curr,0,nChAc);
        }
        iFila.l[curr] = iFila.c;
        (iFila.c)++;
    }
  }
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) {
	nIFilaP = (int)((nActGoal - iFila.c) / nCpu);
	for(n = 0; n < nIFilaP; n++) {
		iFilaP.l[n] = iFila.c + rank * nIFilaP + n;
	}
	iFilaP.c = nIFilaP;
  }
  nAbpM = 0;
  FOR_ABP(n) {
   	pArr = &P2A(chAbp,n,0,nChAb);
    // Count the number of free ABPs
    if (ISMTF(pArr[2])) {
		if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) {
        	nAbpM++;
		}
	}
	else {
		if (pArr[0] < 0 && pArr[1] < 0) {
        	nAbpM++;
		}
	}
  }
  if (nActM % nCpu == 0) { nActM = nActM / nCpu; }
  else { nActM = (int)(nActM / nCpu) + 1; }
  if (nAbpM % nCpu == 0) { nAbpM = nAbpM / nCpu; }
  else { nAbpM = (int)(nAbpM / nCpu) + 1; }

  nActFilaMe = 0;
  nActMe = 0;
  nActCp = 0;
  nAbpMe = 0;
  nAbpCp = 0;
  // If actin is located within the current subdomain, they are listed
  // in iAct and act.id. For example, the actin having the absolute index of 
  // 10 is found for the first time as a particle belonging to the current
  // subdomain, iAct[0] = 10 and act.id[10] = 0. In other words,  
  //   act.id[relative index] = absolute index
  //   iAct[absolute index] = relative index
  // Also, rAct and fixAct are copied to local variables, act.r and act.fix.
  for(m = 0; m < 2; m++) {
	FOR_ACT(n) {
        CONT(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0);
		CS = 0; 
		CS2 = 0;
		ExtractConfigSubroutine(&P2(rAct,n,0), rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) || 
				(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			if (m == 0) { ind = nActMe; }
			else { ind = nActMe + nActCp; }
			V3COPY(&P2(act.r,ind,0), &P2(rAct,n,0));
			V3SET_ALL(&P2(act.fBr,ind,0), 0.);
			act.id[ind] = n;
			act.iF[ind] = iFila.l[n];
			iAct[n] = ind;
			if (m == 0) { 
				nActMe++; 
				if (fixAct[n] > -1) { 
					act.fix[ind] = fixAct[n]; 
				}
                if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
    	        	act.nFA[ind] = (fixAct[n] > -1) ? 1 : 0;
	            }
			}
			else { nActCp++; }
		}	
	}
    CONT(m == 1);
    cnt = 0;
    FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0));
		if (cnt >= nActM * rank && cnt < nActM * (rank + 1)) {
			iAct[n] = nActMe;
			UpdateActinAbpMonomerListSubroutine2(iAct[n], n, 0);
			InsertElement1dArrayWoChk(actM.l, &actM.c, n);
			nActMe++;
		}
		cnt++;
    }
  }
  // Perform the same thing as above for ABPs.
  for(m = 0; m < 2; m++) {
	for(n = 0; n < nAbp; n++) {
   		pArr = &P2A(chAbp,n,0,nChAb);
        CONT((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0));
		CS = 0; 
		CS2 = 0;
		ExtractConfigSubroutine(&P2(rAbp,n,0), rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) || 
					(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			if (m == 0) { ind = nAbpMe; }
			else { ind = nAbpMe + nAbpCp; }
			V3COPY(&P2(abp.r,ind,0), &P2(rAbp,n,0));
			V3SET_ALL(&P2(abp.fBr,ind,0), 0.);
			abp.id[ind] = n;
			iAbp[n] = ind;
			if (m == 0) { 
				nAbpMe++; 
				if (motSA.gTgl != 0) { abp.mId[ind] = abpMotId[n]; }
			}
			else { nAbpCp++; }
		}
	}
    CONT(m == 1);
    cnt = 0;
    FOR_ABP(n) {
    	pArr = &P2A(chAbp,n,0,nChAb);
        CONT(!((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0)));
		if (cnt >= nAbpM * rank && cnt < nAbpM * (rank + 1)) {
			iAbp[n] = nAbpMe;
			UpdateActinAbpMonomerListSubroutine2(iAbp[n], n, pArr[2] + 1);
			if (pArr[2] == 2 && gTglImpMotM != 0) {
				InsertElement1dArrayWoChk(motM.l, &motM.c, n);
			}
			if (pArr[2] != 2 && gTglImpAcpM != 0) {
				InsertElement1dArrayWoChk(acpM.l, &acpM.c, n);
			}
			if ((pArr[2] == 2 && gTglImpMotM == 0) 
					|| (pArr[2] != 2 && gTglImpAcpM == 0)) {
				if (P2(rAbp,n,0) == 0. && P2(rAbp,n,1) == 0. 
						&& P2(rAbp,n,2) == 0.) {
					GenRandPosSubdom(&P2(abp.r,iAbp[n],0));
				}
				else {
					V3COPY(&P2(abp.r,iAbp[n],0), &P2(rAbp,n,0));
				}
			}
			nAbpMe++;
		}
		cnt++;
    }
  }
  if (nActMe > nActC) {
	Printf0("Error: the value of nActC is too small!!\n\n");
	exit(-1);
  }
  if (nAbpMe > nAbpC) {
	Printf0("Error: the value of nAbpC is too small!!\n\n");
	exit(-1);
  }
  nActMeCp = nActMe + nActCp;
  nAbpMeCp = nAbpMe + nAbpCp;
  // Copy chain information from global variables to local variables
  FOR_ACT(n) {
	pArr = &P2A(chAct,n,0,nChAc);
    CONT(iAct[n] < 0 || (pArr[0] < 0 && pArr[1] < 0));
    if (pArr[0] > -1 && pArr[1] < 0 && iAct[n] < nActMe)
    { nActFilaMe++; }
	Copy1dArrayInt(&P2A(act.ch,iAct[n],0,nChAc), pArr, nChAc);
  }
  FOR_ABP(n) {
    pArr = &P2A(chAbp,n,0,nChAb);
    CONT(iAbp[n] < 0);
	CONT((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
			&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
			&& pArr[1] < 0));
	Copy1dArrayInt(&P2A(abp.ch,iAbp[n],0,nChAb), pArr, nChAb);
  }
  if (rheoWay > 0) {
	meaStreParMe.c = 0;
	for(n = 0; n < meaStrePar.c; n++) {
		ind = meaStrePar.l[n];
		if (iAct[ind] > -1 && iAct[ind] < nActMe) {
			meaStreParMe.l[meaStreParMe.c] = ind;
			meaStreParMe.c++;
		}
	}
	appStraParMe.c = 0;
	for(n = 0; n < appStraPar.c; n++) {
		ind = appStraPar.l[n];
		if (iAct[ind] > -1 && iAct[ind] < nActMe) {
			appStraParMe.l[appStraParMe.c] = ind;
			appStraParMe.c++;
		}
	}
  }
  // Free unnecessary global arrays.
  free(rAct);
  free(rAbp);
  free(chAct);
  free(chAbp);
  free(fixAct);
  if (motSA.gTgl != 0) { free(abpMotId); }
  if (recTraj.gTglCho == 0) {
	free(recTraj.act.l); 
	free(recTraj.abp.l); 
  }
  if (rheoWay > 0) {  
	free(appStraPar.l);
	free(meaStrePar.l); 
  }
  free(iFila.l);
  UpdateNeighborList();
  UpdateChainList();

}

// Load information about network configuration from "Config".
void LoadConfig(int mode) {
  int n, k, l, ind, tempInt, tempInt2[2], tglTraj, tglFix;
  char ch[200];
  FILE *fIn; 

  nMot = 0;
  // mode 0 is for loading data to resume simulation in case of errors
  // mode 1 is for the initial loading of configuration file
  if ((fIn = fopen("Config", "r")) == NULL) {
	Printf0("File doesn't exist: Config\n");
	exit(-1);
  }
  // nAct, nAbp, are dimDom[] are loaded at LoadInitParameter(), not here. 
  while (strcmp(ch, "## Position for actin ##\n")) {
	fgets(ch, 200, fIn); 
  }

  // Positions of actin and ABP
  FOR_ACT(n) {
    fscanf(fIn, "%d\t%lf\t%lf\t%lf", &tempInt, &P2(rAct,n,0), &P2(rAct,n,1),
		&P2(rAct,n,2));
  }

  fscanf(fIn, "\n## Position for ABP ##\n");
  FOR_ABP(n) {
	fscanf(fIn, "%d\t%lf\t%lf\t%lf", &tempInt, &P2(rAbp,n,0), &P2(rAbp,n,1), 
		&P2(rAbp,n,2));
  }

  // Chain information of actin and ABP
  memset(chAct, -1, sizeof(int) * nAct * nChAc);
  fscanf(fIn, "\n## Chain for actin ##\n");
  FOR_ACT(n) {
	fscanf(fIn, "%d", &tempInt);
	for(k = 0; k < 2; k++) {
		fscanf(fIn, "%d", &P2A(chAct,n,k,nChAc));
	}
	for(k = 0; k < confNChAcX; k++) {
		for(l = 0; l < confNChAcY; l++) {
			ind = 2 + k * nChAcY * (nChAcX / confNChAcX) + l;
			fscanf(fIn, "%d", &P2A(chAct,n,ind,nChAc));
		}
	}
  }

  fscanf(fIn, "\n## Chain for ABP ##\n");
  FOR_ABP(n) {
	fscanf(fIn, "%d", &tempInt);
	for(k = 0; k < nChAb; k++) {
		fscanf(fIn, "%d", &P2A(chAbp,n,k,nChAb));
	}
	if (P2A(chAbp,n,2,nChAb) == 2) { nMot++; }
  }

  fscanf(fIn, "\n## Position for membrane ##\n");
  FOR_MB(n) {
	fscanf(fIn, "%d", &tempInt);
	fscanf(fIn, "%lf\t%lf\t%lf", &P2(rMb,n,0), &P2(rMb,n,1), 
			&P2(rMb,n,2));
  }
  fscanf(fIn, "\n## Chain for membrane ##\n");
  FOR_MB(n) {
	fscanf(fIn, "%d", &tempInt);
	for(k = 0; k < nChMb; k++) {
		fscanf(fIn, "%d", &P2A(chMb,n,k,nChMb));
		if (k >= nChMb - nMbAct && P2A(chMb,n,k,nChMb) > -1) {
			actRebMb[P2A(chMb,n,k,nChMb)] = n;
		}
	}
  }
  // Related to rheological measuremenets
  tglTraj = (recTraj.gTglCho == 0) ? 1 : 0;
  tglFix = (rheoWay > 0 && gTglLoadNetDataFix != 0) ? 1 : 0;

  while(!(ch[0] == 's' && ch[1] == 't')) {
	fgets(ch, 200, fIn);
  }
  fgets(ch, 200, fIn);
  fscanf(fIn, "meaStrePar = %d\n", &meaStrePar.c);
  for(n = 0; n < meaStrePar.c; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglFix != 0) { meaStrePar.l[n] = tempInt; }
  }
  fscanf(fIn, "\nappStraPar = %d\n", &appStraPar.c);
  for(n = 0; n < appStraPar.c; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglFix != 0) { appStraPar.l[n] = tempInt; }
  }
  fscanf(fIn, "\nactTraj = %d\n", &recTraj.act.c);
  if (tglTraj != 0) { recTraj.nActL = recTraj.act.c; }
  for (n = 0 ; n < recTraj.act.c ; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglTraj != 0) { recTraj.act.l[n] = tempInt; }
  }
  fscanf(fIn, "\nabpTraj = %d\n", &recTraj.abp.c);
  if (tglTraj != 0) { recTraj.nAbpL = recTraj.abp.c; }
  for (n = 0 ; n < recTraj.abp.c ; n++) {
	fscanf(fIn, "%d\n", &tempInt);
	if (tglTraj != 0) { recTraj.abp.l[n] = tempInt; }
  }
  fscanf(fIn, "\nfixAct = %d\n", &tempInt);
  for(n = 0; n < tempInt; n++) {
	fscanf(fIn, "%d\t%d\n", &tempInt2[0], &tempInt2[1]);
	if (gTglLoadNetDataFix != 0) {
		fixAct[tempInt2[0]] = tempInt2[1];
	}
  }
  fclose(fIn);
}

// Subroutine for LoadInitParameter(), which judges "yes" or "no" from 
// the file, "condition".
void CheckAnswerYesNo(char *tag, const char *which, int *tgl) {
  if (strcmp(tag, "yes") == 0) { *tgl = 1; }
  else if (strcmp(tag, "no") == 0) { *tgl = 0; } 
  else { 
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

/*------------------------ Related to network data ---------------------------*/

/*----------------------- Loading initial parameters -------------------------*/

void LoadInitParameterSubroutine2(char *str, int ind, int *tgl1, int *tgl2) {
  *tgl1 = (ind > 0) ? 1 : 0;
  *tgl2 = (ind == 2) ? 1 : 0;
  if (ind < 0 || ind > 2) {
	Printf0("Error: choice of %s should be either of 0, 1, or 2.\n\n", str);
  	MPI_Barrier(MPI_COMM_WORLD);
	exit(-1);
  }
}

void LoadInitParameterSubroutine(FILE *fIn, const char *which, int *tgl) {
  char tag[80];

  fscanf(fIn, "%s", tag);
  if (tag[0] == 'y' && tag[1] == 'e' && tag[2] == 's') 
  { *tgl = 1; }
  else if (tag[0] == 'n' && tag[1] == 'o') 
  { *tgl = 0; }
  else {
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

// Load parameters from "Config" and "condition" at the very beginning.
// Those parameters are needed to define arrays, and so on.
void LoadInitParameter(void) {
  int n, k, tempInt[4], CS;
  double tempDbl[3];
  char tag[200], direc[4]= "xyz";
  FILE *fIn;

  nChAb = 5;
  MALLOC(bndUnb.facK0,double,NDIM*2);
  MALLOC(bndUnb.facX,double,NDIM*2);
  MALLOC(bndReb.facK,double,NDIM*2);
  MALLOC(bnd.drag,double,NDIM*2);
  MALLOC(bnd.stfSpr,double,NDIM*2);
  MALLOC(abpF.a90,Force,NK_ABP);
  MALLOC(abpF.cr,Force,NK_ABP);
  MALLOC(abpF.bend,Force,NK_ABP);
  MALLOC(abpF.spr,Force,NK_ABP);

  // From "condition", several types of informatio are loaded: MPI method,
  // rheological methods, dynamic behaviors of ABPs
  if ((fIn = fopen("condition", "r")) == NULL) {
	Printf0("File doesn't exist: condition\n");
	exit(-1); 
  }
  fgets(tag, 200, fIn);

  mpiMethod = 0.;
  updSubdSize.prdR = -1;
  nActPerSeg = 20;
  nChAcX = 20;
  nChAcY = 10;

  // Network condition
  fscanf(fIn, "======================= Parameters for network condition "
		"=======================\n");
  fscanf(fIn, "Load network data from Config(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "load network data", &gTglLoadNetData);
  fscanf(fIn, "If no, width of domain(x, y, z in um) = %lf, %lf, %lf\n", 
		&tempDbl[0], &tempDbl[1], &tempDbl[2]);
  if (gTglLoadNetData == 0) {
	FOR_NDIM(n) { dimDom[n] = L_UM2S(tempDbl[n]); }
  }
  fscanf(fIn, "Periodic boundary condition"
		"(yes/no in x, y, z) = ");
  FOR_NDIM(n) {
	LoadInitParameterSubroutine(fIn, "PBC", &pbc[n]);
  }
  fscanf(fIn, "\n");
  fscanf(fIn, "If 2D network, specifiy the normal direction(no, x, y, or z) "
		"= %s\n", tag);
  if (strcmp(tag, "x") == 0) { dir2D = 0; }
  else if (strcmp(tag, "y") == 0) { dir2D = 1; } 
  else if (strcmp(tag, "z") == 0) { dir2D = 2; } 
  else if (strcmp(tag, "no") == 0) { dir2D = -1; } 
  else {
	Printf0("Error: Normal direction specifying 2D geometry should be "
			"either no, x, y, or z.\n\n"); 
	exit(-1); 
  }
  if (dir2D > -1) {
	if (pbc[dir2D] == 1) {
		Printf0("Error: periodic boundary condition should not exist "
				"in %c direction with the 2D geometry.\n\n", direc[dir2D]); 
		exit(-1); 
	}
  }

  gTglLoadNetDataFix = 0;

  fscanf(fIn, "Actin concentration(in uM or given) = %s\n", tag);
  if (strcmp(tag, "given") == 0) { 
	if (gTglLoadNetData == 0) {
		Printf0("Error: no given actin concentration without "
				"loaded network data.\n\n"); 
		exit(-1); 
	}
	cAct = -1.; 
  }
  else { cAct = atof(tag); }
  // R values of ACPC, ACPB, and motor
  fscanf(fIn, "ACPC density(R value or given) = %s\n", tag);
  if (strcmp(tag, "given") == 0) { RAbp[0] = -1.; }
  else { RAbp[0] = atof(tag); }
  fscanf(fIn, "ACPB density(R value or given) = %s\n", tag);
  if (strcmp(tag, "given") == 0) { RAbp[1] = -1.; }
  else { RAbp[1] = atof(tag); }
  fscanf(fIn, "Motor density(R value or given) = %s\n", tag);
  if (strcmp(tag, "given") == 0) { RAbp[2] = -1.; }
  else { RAbp[2] = atof(tag); }
  if (gTglLoadNetData == 0 && (RAbp[0] == -1. || RAbp[1] == -1. 
		|| RAbp[2] == -1.)) {
	Printf0("Error: no given ABP concentration without "
			"loaded network data.\n\n"); 
	exit(-1); 
  }

  fscanf(fIn, "Duration of network formation(s) = %lf\n", &netForm.durR);
  fscanf(fIn, "Duration of simulation(s) = %lf\n", &rheo.durR);

  rheoWay = -1;
  motActiv.durR = 0.;
  recTraj.gTglCho = 0;
  recTraj.porActL = 0.;
  recTraj.porAbpL = 0.;
  recTraj.prdR = 1.;

  bulkRheoWay = 1;
  bulkRheoType = 0;
  signStr = 1;
  dirStr = 0;
  dirNoPBC = 2;
  pres.mag = 0.;
  pres.rate = 1e5;
  sinuStr.amp = 1e-6;
  sinuStr.prdR = 1.;
  stra.lim[0] = -10.;
  stra.lim[1] = 10.;
  recStre.prdR = 1e-4;
  gTglBead = 0;
  bead.n = 0;
  bead.rad2 = 0.5;
  bead.facStfRep = 1.;
  bead.thkRep = 0.05;
  recBeadLoc.prdR = 0.01;
  beadBind.gTgl = 0;
  bead.rheoWay = 1;
  bead.rheoDir = 0;
  bead.rheoSig = 0;
  bead.rheoPreMag = 0.;
  bead.rheoPreRate = 10.;
  bead.rheoAmp = 10.;
  bead.rheoPrdR = 1.;

  // Toggle dynamic behaviors of actin and ABP
  fscanf(fIn, "==================== Toggle the dynamic behaviors of actin "
		"=====================\n");
  fscanf(fIn, "Allow the thermal fluctuation of actin filaments(yes/no) "
		"= %s\n", tag);
  CheckAnswerYesNo(tag, "actin thermal fluctuation", &gTglActTherm);

  gTglActDynNF = 1;
  gTglActDynPres = 1;

  fscanf(fIn, "Allow nucleation of actins(0: no, 1: stochastic, 2: constant "
		"# of filaments) = %d\n", &tempInt[0]);
  LoadInitParameterSubroutine2("actin nucleation", tempInt[0], 
			&actNuc.gTgl, &actNuc.gTglFN);
  if (actNuc.gTglFN != 0 && gTglLoadNetData == 0) {
	Printf0("Error: the number of actin filaments cannot be maintained because"
			"network data are not load from a file!");
	exit(-1); 
  }

  actBst.gTgl = 0;
  actCap.gTgl = 0;
  actUnc.gTgl = 0;
  actSev.gTgl = 0;
  actSev.gTglCap[0] = 0;
  actSev.gTglCap[1] = 0;
  actSev.gTglBst[0] = 0;
  actSev.gTglBst[1] = 0;
  actAnn.gTgl = 0;
  actDgd.gTgl = 0;

  fscanf(fIn, "Allow assembly of actins(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin assembly", &actAss.gTgl);
  fscanf(fIn, "Allow disassembly of actins(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin disassembly", &actDis.gTgl);

  gTglAcpDynNF = 1;
  gTglAcpDynPres = 1;

  fscanf(fIn, "===================== Toggle the dynamic behaviors of ACP "
		"======================\n");
  fscanf(fIn, "Allow the thermal fluctuation of ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "ACP thermal fluctuation", &gTglAcpTherm);
  fscanf(fIn, "Allow unbinding of inactive ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "inactive ACP unbinding", &acpInaUnb.gTgl);
  fscanf(fIn, "Allow binding of monomeric ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "monomeric ACP binding", &acpMoBind.gTgl);
  fscanf(fIn, "Allow unbinding of active ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "active ACP unbinding", &acpUnb.gTgl);
  fscanf(fIn, "Allow binding of inactive ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "inactive ACP binding", &acpReb.gTgl);

  acpReb.gTglCrsAng = 0;
  gTglImpAcpM = 1;

  fscanf(fIn, "==================== Toggle the dynamic behaviors of motor "
		"=====================\n");
  fscanf(fIn, "Allow the thermal fluctuation of motors(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "motor thermal fluctuation", &gTglMotTherm);
  fscanf(fIn, "Allow walking of motors(0: no, 1: w/o slide-off, "
		"2: w/ slide-off) = %d\n", &tempInt[0]);
  LoadInitParameterSubroutine2("motor walking", tempInt[0], 
		&motWalk.gTgl, &gTglMotWalkSld);

  gTglMotUnbRebNF = 1;
  gTglMotUnbRebPres = 1;
  gTglMotWalkNF = 0;
  gTglMotWalkPres = 1;
  motSA.gTgl = 1;
  motSA.gTglConSiz = 1;
  motReb.gTglCrsAng = 0;
  gTglImpMotM = 1;

  fscanf(fIn, "Allow unbinding of motors(yes/no) = %s\n", tag);
  if (motSA.gTgl != 0) {
	CheckAnswerYesNo(tag, "motor unbinding", &motInaUnb.gTgl);
	motUnb.gTgl = motInaUnb.gTgl;
  }
  fscanf(fIn, "Allow binding of motors(yes/no) = %s\n", tag);
  if (motSA.gTgl != 0) {
	CheckAnswerYesNo(tag, "motor binding", &motMoBind.gTgl);
	motReb.gTgl = motMoBind.gTgl;
  }

  motReb.gTglCrsAng = 1;
  motReb.gTglOppDir = 0;
  motSA.to.gTgl = 1;
  motSA.to.mode = 0;
  updMono.prdR = 0.1;

  // Parameters for dynamic behaviors of actin and ABP
  fscanf(fIn, "================ Parameters for the dynamic behaviors of "
		"actin =================\n");
  fscanf(fIn, "k for actin nucleation(1/uM s) = %lf\n", &actNuc.k);
  fscanf(fIn, "k for actin assembly at barbed and pointed ends(1/uM s) "
		"= %lf, %lf\n", &actAss.k[0], &actAss.k[1]);
  fscanf(fIn, "k for actin disassembly at barbed and pointed ends(1/s) "
		"= %lf, %lf\n", &actDis.k[0], &actDis.k[1]);
  fscanf(fIn, "Factor for varying disassembly rate of actin with ACPs or "
		"motors = %lf\n", &actDis.facKWA);

  actBst.k[0] = 0.1;
  actBst.k[1] = 0.1;
  actBst.facKWA = 1.;
  actSev.k = 1e-40;
  actSev.facX = 1.6;
  actSev.critF = 500;
  actSev.facKWA = 1.;
  actAnn.k = 1.;
  actAnn.ang = 10.;
  actCap.k[0] = 1.;
  actCap.k[1] = 1.;
  actUnc.k[0] = 1.;
  actUnc.k[1] = 1.;
  actDgd.k = 10.;
  actDgd.x1 = 0;
  actDgd.x2 = 0;
  actDgd.dist = 150;

  fscanf(fIn, "================= Parameters for the dynamic behaviors of "
		"ACP ==================\n");
  fscanf(fIn, "k for ACP binding(k=i*k_0) = %lf\n", &acpReb.facK);
  acpMoBind.facK = acpReb.facK;
  fscanf(fIn, "k0 and x for ACP unbinding(k0=i1*k0_0, "
		"x=i2*x_0) = %lf, %lf\n", &acpUnb.facK0, &acpUnb.facX);
  fscanf(fIn, "================ Parameters for the dynamic behaviors of "
		"motor =================\n");
  fscanf(fIn, "k for motor binding(k=i*k_0) = %lf\n", 
		&motReb.facK);
  motMoBind.facK = motReb.facK;
  fscanf(fIn, "k0 for motor unbinding and walking(k0=i*k0_0) = %lf\n", 
		&motUnb.facK0);
  motWalk.facK0 = motUnb.facK0;
  fscanf(fIn, "Number of heads which each motor arm represents = %d\n", 
		&motMC.nHead);
  fscanf(fIn, "Transition rates between mechanochemical states in each head"
		"(k01, k10, k12, k21, k20 in 1/s) = %lf, %lf, %lf, %lf, %lf\n", 
		&motMC.k01, &motMC.k10, &motMC.k12, &motMC.k21, &motMC.k20);
  fscanf(fIn, "Average number of motors per each self-assembled structure = "
		"%d\n", &motSA.nMotPerTF);
  fscanf(fIn, "Portion of inactive motor arms in self-assembled structures "
		"= %lf\n", &motSA.porInac);

  updMotN.prdR = 1.;

  fscanf(fIn, "k for motor assembly(in 1/s) = %lf\n", &motSA.kAss);
  fscanf(fIn, "k for motor turnover(in 1/s) = %lf\n", &motSA.to.k);

  motSA.to.k0 = 1e-3;
  motSA.to.x = -1e-10;

  fscanf(fIn, "======== Adjustment for rates of the dynamic behaviors of "
		"actin and ABP ========\n");
  fscanf(fIn, "During network formation(i*k) = %lf\n", &netForm.facK);

  pres.facK = 1.;
  motActiv.facK = 1.;

  fscanf(fIn, "=============== Parameters for the mechanical behaviors of "
		"actin ===============\n");
  fscanf(fIn, "Strength of repulsive forces between actin filaments(Kr=i*Kr_0)"
		" = %lf\n", &actF.rep.facStf);
  fscanf(fIn, "Bending stiffness of actin filaments(Kb=i*Kb_0) = %lf\n", 
		&actF.bend.facStf);
  fscanf(fIn, "Extensional stiffness of actin filaments(Ks=i*Ks_0) = %lf\n", 
		&actF.spr.facStf);

  fscanf(fIn, "================ Parameters for the mechanical behaviors of "
		"ABP ================\n");

  abpF.rep.facStf = 0.;

  fscanf(fIn, "Bending stiffness which maintains an assigned angle between two "
		"arms of ACPC, ACPB, and motor(Kb=i*Kb_0) = %lf, %lf, %lf\n", 
		&abpF.bend[0].facStf, &abpF.bend[1].facStf, &abpF.bend[2].facStf);
  fscanf(fIn, "Bending stiffness which maintains 90 deg between the axis of "
		"a filament and the arm of ACPC, ACPB, and motor(Kb=i*Kb_0) = %lf, "
		"%lf, %lf\n", &abpF.a90[0].facStf, &abpF.a90[1].facStf, 
		&abpF.a90[2].facStf);
  fscanf(fIn, "Extensional stiffness of ACPC, ACPB, and motor(Ks=i*Ks_0) = "
		"%lf, %lf, %lf\n", &abpF.spr[0].facStf, &abpF.spr[1].facStf, 
		&abpF.spr[2].facStf);

  bnd.radRnd = -1;
  bndMv.gTgl = 0;
  for(k = 0; k < 6; k++) {
	bndMv.stf[k] = 100.;
	bndMv.thk[k] = 5.;
  }
  bndMv.stfUnit = 0;
  bndVol.gTgl = 0;
  bndVol.facStf = 0.;

  bndUnb.gTgl = 0;
  for(k = 0; k < 6; k++) {
	bndUnb.facK0[k] = 1.;
	bndUnb.facX[k] = 1.;
	bndReb.facK[k] = 1.;
	bnd.stfSpr[k] = 1.;
  }
  bndReb.gTgl = 0;
  bndReb.gTglPa = 0;
  bndReb.depR = 50.;
  bndReb.gTglSpr = 0;
  recBndLoc.prdR = -1;
  bnd.gTglActMv = 0;
  for(k = 0; k < 6; k++) {
	bnd.drag[k] = 100.;
  }
  bndMat.gTgl = 0;
  bndMat.gTglPa = 0;
  bndMat.facK0 = 1.;
  bndMat.facX = 1.;
  bndMat.maxNFA = 20;
  recBndUnbReb.prdR = -1;
  recBndActMat.prdR = -1;
  recBndTracF.prdR = -1;
  gTglMb = 0;
  gTglNuc = 0;
  gTglLoadMbNucData = 0;
  dimMbNuc = 3;
  dirNormMbNuc = 2;
  nObjMbNuc = 1;
  radMb = 3.0;
  memb.len = 0.35;
  memb.thk = 0.1;
  memb.dragR = 1.0;
  memb.bend.facStf = 1.0;
  memb.spr.facStf = 1.0;
  memb.spr2.facStf = 1.0;
  mbVol.gTgl = 1;
  memb.vol.facStf = 1.;
  mbAre.gTgl = 1;
  memb.area.facStf = 1.0;
  gTglMbTherm = 1;
  mbPro.gTgl = 0;
  mbPro.facK = 1.0;
  mbPro.gTglCT = 0;
  mbPro.rCT[0] = 0.5;
  mbPro.rCT[1] = mbPro.rCT[2] = 0.;
  mbPro.f = 100.;
  mbPro.durR = 5.;
  mbFix.gTgl = 0;
  mbFix.facK = 0.01;
  mbDef.gTgl = 0;
  mbDef.facK0 = 1.;
  mbDef.facX = 1.;
  sideMb = 0;
  thkMbActNuc = 0.1;
  gTglMbCont = 0;
  mbReb.gTgl = 0;
  mbSld.gTgl = 0;
  mbReb.gTglPa = 0;
  mbReb.por = 1.0;
  nMbAct= 5;
  mbReb.facK = 1.0;
  mbUnb.gTgl = 0;
  mbUnb.facK0 = 1.;
  mbUnb.facX = 1.;
  mbMat.gTgl = 0;
  mbMat.gTglPa = 0;
  mbMat.facK0 = 1.;
  mbMat.facX = 1.;
  mbMat.maxNFA = 20;
  mbSld.act.gTgl = 0;
  mbSld.abp.gTgl = 0;
  mbSld.act.por = 0.1;
  mbSld.abp.por = 0.1;
  mbSld.critDist = 300.;
  mbSld.eqDist = 100;
  mbSld.stf = 0.0003;
  mbSld.drag = 10.;
  radNuc = 2.;
  memb.nucLen = 0.35;
  memb.nucThk = 0.2;
  memb.nucDragR = 1.;
  memb.nucBend.facStf = 10.;
  memb.nucSpr.facStf = 1.;
  memb.nucSpr2.facStf = 1.;
  nucVol.gTgl = 1;
  memb.nucVol.facStf = 6.;
  nucAre.gTgl = 0;
  memb.nucArea.facStf = 1.;
  gTglNucTherm = 0;
  nucReb.gTgl = 0; 
  nucReb.gTglPa = 0;
  nucReb.por = 1.;
  nNucAct = 5;
  nucReb.facK = 1.;
  nucUnb.gTgl = 0;
  nucUnb.facK0 = 1.;
  nucUnb.facX = 1.;

  fscanf(fIn, "======================== Parameters for data recording ======"
		"===================\n");
  fscanf(fIn, "Period of recording Output and Progress(s) = %lf\n", 
		&recProg.prdR);

  fscanf(fIn, "----------------------------- Structural information "
		"---------------------------\n");
  fscanf(fIn, "Period of recording Config(in s, -1 for deactivation) = %lf\n",
		&recConf.prdR);
  fscanf(fIn, "Period of recording structural information for visualization "
		"via VMD(in s, -1 for deactivation) = %lf\n", &recConfVmd.prdR);
  fscanf(fIn, "Number of VMD files(multiple/single) = %s\n", tag);
  if (strcmp(tag, "multiple") == 0) { recConfVmd.mode = 0; }
  else if (strcmp(tag, "single") == 0) { recConfVmd.mode = 1; } 
  else {
	Printf0("Error: the number of VMD files should be multiple or single\n");
	exit(-1); 
  }


  recConfMlb.prdR = -1;
  recMotSize.prdR = -1;
  recMotPos.prdR = -1;
  recPoreSize.prdR = -1;
  recConn.prdR = -1;
  recPerc.prdR = -1;
  dirRecPerc = 2;
  findSupp.prdR = -1;
  porFindSupp = 0.2;
  kindFindSupp = 0;

  fscanf(fIn, "Show the boundaries of a network drawn by VMD(yes/no) = %s\n", 
		tag);
  CheckAnswerYesNo(tag, "show boundaries in VMD", &recConfVmd.gTglBnd);
  fscanf(fIn, "Record information for coloring a network drawn by "
		"VMD(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "information for coloring VMD", &recConfVmd.gTglInfo);
  fscanf(fIn, "Lower and upper limits of forces for coloring(in pN) = %lf, "
		"%lf\n", &recConfVmd.minF, &recConfVmd.maxF);
  fscanf(fIn, "Period of recording length of actin filaments(in s, -1 for "
		"deactivation) = %lf\n", &recFilaL.prdR);
  fscanf(fIn, "Period of recording distances between active ABPs(in s, -1 for"
		" deactivation) = %lf\n", &recCrsDist.prdR);

  fscanf(fIn, "--------------------------- Force, stress, and energy "
		"--------------------------\n");
  fscanf(fIn, "Period of recording longitudinal forces of ABPs(in s, -1 for "
		"deactivation) = %lf\n", &recLongF.prdR);
  fscanf(fIn, "Period of recording mechanical energy(in s, -1 for "
		"deactivation) = %lf\n", &recE.prdR);
  fscanf(fIn, "Period of recording internal elastic/viscous stresses"
		"(in s, -1 for deactivation) and the number of measurements("
		"in x, y, z) = %lf, %d, %d, %d\n", &recSecStre.prdR, 
		&nSecStreDiv[0], &nSecStreDiv[1], &nSecStreDiv[2]);

  recSecStre.mode = 0;
  recActSev.prdR = -1;
  recAbpUnb.prdR = -1;
  recAbpBind.prdR = -1;
  recAbpTurn.prdR = -1;
  recAbpDyn.prdR = -1;
  recMbCen.prdR = -1;
  recMbDim.prdR = -1;

  fscanf(fIn, "---------------------------------- Miscellany "
		"----------------------------------\n");
  fscanf(fIn, "Period of recording information in unit of individual "
		"elements, filament segments, and filaments(in s, -1 for deactivation)"
		" = %lf\n", &recInfo.prdR);
  fclose(fIn);

  if (cAct >= 0. && actAss.gTgl == 0 && actDis.gTgl == 0 
		&& actBst.gTgl == 0) { 
	cAct = -1.; 
  }
  if (gTglLoadNetData == 0 && actAss.gTgl == 0 && actNuc.gTgl == 0) {
	Printf0("Error: actin assembly and nucleation should be allowed "
			"at least.\n\n"); 
	exit(-1); 
  }
  if (rheoWay > 0 && bndMv.gTgl != 0) {
	Printf0("Warning: bulk rheology cannot be used with moving boundaries. "
			"The moving boundaries will be deactivated!!\n\n");
	bndMv.gTgl = 0;
  }
  if (rheoWay == 0 || rheoWay == 2) {
	if ((bndUnb.gTgl != 0 || bndReb.gTgl != 0 || bnd.gTglActMv != 0) 
			&& pbc[0] != 0 && pbc[1] != 0 && pbc[2] != 0) {
		Printf0("Warning: periodic boundary conditions exist in all directions,"
				" so the unbinding/binding of actins on boundaries will be"
				" ignored!\n\n");
		bndUnb.gTgl = bndReb.gTgl = bnd.gTglActMv = 0;
	}
  }

  V3SET_ALL(nAbpDet, 0);
  if (gTglLoadNetData != 0) {
	// From "Config", the values of nAct, nAbp, and dimDom[] are loaded.
	if ((fIn = fopen("Config", "r")) == NULL) {
	    Printf0("File doesn't exist: Config\n");
	    exit(-1);
	}
	fscanf(fIn, "nAct : %d\n", &nAct);
	fscanf(fIn, "nAbp : %d\n", &nAbp);
	fscanf(fIn, "nMb : %d\n", &nMb);
	fscanf(fIn, "dimDom : %lf, %lf, %lf\n", &dimDom[0], &dimDom[1], &dimDom[2]);

	fscanf(fIn, "nActPerSeg, nChAcX, nChAcY : %d, %d, %d\n", &tempInt[0], 
			&confNChAcX, &confNChAcY);
	if (tempInt[0] != nActPerSeg){
	    Printf0("Error: the length of actin cylindrical segments "
				"doesn't match between 'condition' and 'Config': "
				"%g nm in 'condition' vs %g nm in 'Config'\n", 
				nActPerSeg * 7., tempInt[0] * 7.);
	    exit(-1);
	}
	if (confNChAcX > nChAcX) {
		Printf0("Error: the number of binding sites in a longitudinal direction"
				" on each actin segment in 'Config' is greater than that in"
				" 'condition': %d in 'Config' vs %d in 'condition'\n", 
				confNChAcX, nChAcX);
	    exit(-1);
	}
	else if (confNChAcX < nChAcX) {
		if (nChAcX % confNChAcX == 0) {
			Printf0("Warning: the number of binding sites in a longitudinal"
					" direction on each actin segment in 'Config' is smaller"
					" than that in 'condition': %d in 'Config' vs %d in"
					" 'condition'\n\n", confNChAcX, nChAcX);
		}
		else {
			Printf0("Error: the number of binding sites in a longitudinal"
					" direction on each actin segment in 'Config' is smaller"
					" than that in 'condition'. In this case, the latter"
					" should be a multiple of the former: %d in 'Config'"
					" vs %d in 'condition'\n", confNChAcX, nChAcX);
		    exit(-1);
		}
	}
	if (confNChAcY > nChAcY) {
		Printf0("Error: the number of binding sites in a transverse direction"
				" on each actin segment in 'Config' is greater than that in"
				" 'condition': %d in 'Config' vs %d in 'condition'\n", 
				confNChAcY, nChAcY);
	    exit(-1);
	}
	else if (confNChAcY < nChAcY) {
		Printf0("Warning: the number of binding sites in a transverse"
				" direction on each actin segment in 'Config' is smaller"
				" than that in 'condition': %d in 'Config' vs %d in"
				" 'condition'\n\n", confNChAcY, nChAcY);
	}

	while (strcmp(tag, "## Chain for actin ##\n")) {
		fgets(tag, 200, fIn); 
	}
	FOR_ACT(n) {
		for(k = 0; k < 3 + confNChAcX * confNChAcY; k++) { 
			fscanf(fIn, "%d", &tempInt[0]);
		}
	}
	while (strcmp(tag, "## Chain for ABP ##\n")) {
		fgets(tag, 200, fIn); 
	}
	FOR_ABP(n) {
		for(k = 0; k < nChAb + 1; k++) { 
			fscanf(fIn, "%d", &tempInt[0]);
			if (k == 3) { nAbpDet[tempInt[0]]++; }
		}
	}
	while(!(tag[0] == 'r' && tag[1] == 'G')) {
		fgets(tag, 200, fIn); 
	}
	fscanf(fIn, "pbc = %d, %d, %d\n", &confPbc[0], &confPbc[1], &confPbc[2]);
  }
  else {
	nAct = 0;
	nAbp = 0;
  }
  CS = 1;
  if (bndUnb.gTgl == 0) {
	free(bndUnb.facK0);
	free(bndUnb.facX);
  }
  if (bndReb.gTgl == 0) {
	free(bndReb.facK);
  }
  if (bnd.gTglActMv == 0) {
	free(bnd.drag);
  }
  gTglActCapAll = (actCap.gTgl != 0 || actUnc.gTgl != 0
		|| (actSev.gTgl != 0 && gTglActSevCap != 0)) ? 1 : 0;
}

/*----------------------- Loading initial parameters -------------------------*/

/*----------------------- Adding and deleting elements -----------------------*/

void AddFreeActinAbpSubroutine(int nMall, int *nAddMe, int *nAddMeSum) {
  int k, mSum, *mAll;
  double volSum, *volAll;

  volMe = HowManyInOutBound();
  MALLOC(volAll, double, nCpu);
  MALLOC(mAll, int, nCpu);
  MPI_Gather(&volMe, 1, MPI_DOUBLE, volAll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	volSum = SumArrDbl(volAll, nCpu);
	for(k = 0; k < nCpu; k++) {
		mAll[k] = (int)(volAll[k] / volSum * nMall);
	}
	mSum = SumArrInt(mAll, nCpu);
	mAll[nCpu - 1] += nMall - mSum;
  }
  MPI_Bcast(mAll, nCpu, MPI_INT, 0, MPI_COMM_WORLD);
  *nAddMe = mAll[rank];
  *nAddMeSum = SumArrInt(mAll, rank);

  free(volAll);
  free(mAll);
}

void AddFreeAbp(void) {
  int n, kind, CS, sum, nAbpAddSumMe, nAbpMallSum, tempInt, tempInt2, cnt;
  int nAbpAddMe[NK_ABP], sum2[NK_ABP], *nAbpAddSumAll;
  int diff, ind, *nNucAll, *cntMotAll;

  int cntMot = 0;
  FOR_ABPME(n) {
	if (K_ABP(n) == 2) { 
		cntMot++; 
	}
  }

  Printf0("%d ACPC, %d ACPB, %d motor are added.\n",
		nAbpMall[0], nAbpMall[1], nAbpMall[2]);
  nMot += nAbpMall[2];
  nAbpMallSum = V3SUM(nAbpMall);

  for(n = 0; n < NK_ABP; n++) {
	CONT(motSA.gTgl != 0 && n == 2);
	nAbpAddMe[n] = nAbpMall[n] / nCpu
			+ (rank < (nAbpMall[n] % nCpu) ? 1 : 0);
  }
  if (motSA.gTgl != 0) {
	MALLOC(nNucAll, int, nCpu * 2);
	MALLOC(cntMotAll, int, nCpu);
	MPI_Gather(&cntMot, 1, MPI_INT, cntMotAll, 1, MPI_INT, 0, 
			MPI_COMM_WORLD);

	if (rank == 0) {
		for(n = 0; n < nCpu; n++) {
			P2A(nNucAll,n,0,2) = motSA.nNucMe / nCpu;
		}
		if (motSA.nNucMe / nCpu > 0) {
			for(n = 0; n < motSA.nNucMe % nCpu; n++) {
				while(1) {
					ind = GenRandIntIndex(nCpu);
					BREAK(P2A(nNucAll,ind,0,2) == motSA.nNucMe / nCpu);
				}
				P2A(nNucAll,ind,0,2)++;
			}
			for(n = 0; n < nCpu; n++) {
				P2A(nNucAll,n,1,2) = P2A(nNucAll,n,0,2) 
						* (motSA.nMotPerSide * 2) - cntMotAll[n];
			}
		}
		else {
			cnt = nAbpDet[2];
			for(n = 0; n < motSA.nNucMe; n++) {
				while(1) {
					ind = GenRandIntIndex(nCpu);
					CONT(cnt > 0 && cntMotAll[ind] == 0); 
					BREAK(P2A(nNucAll,ind,0,2) == 0);
				}
				P2A(nNucAll,ind,0,2)++;
				cnt--;
			}
			for(n = 0; n < nCpu; n++) {
				P2A(nNucAll,n,1,2) = P2A(nNucAll,n,0,2) 
						* (motSA.nMotPerSide * 2) - cntMotAll[n];
			}
		}
	}
	MPI_Bcast(nNucAll, nCpu * 2, MPI_INT, 0, MPI_COMM_WORLD);
	motSA.nNucMe = P2A(nNucAll,rank,0,2);
	nAbpAddMe[2] = P2A(nNucAll,rank,1,2);
	free(nNucAll);
	free(cntMotAll);
  }
  nAbpAddSumMe = SumArrInt(nAbpAddMe, NK_ABP);
  MALLOC(nAbpAddSumAll, int, nCpu);
  MPI_Allgather(&nAbpAddSumMe, 1, MPI_INT, nAbpAddSumAll, 1, MPI_INT,
		MPI_COMM_WORLD);
  sum = SumArrInt(nAbpAddSumAll, rank);
  free(nAbpAddSumAll);

  for(n = 0; n < NK_ABP; n++) {
    nAbpDet[n] += nAbpMall[n];
    nAbpMall[n] = 0;
  }
  nAcpMme += nAbpAddMe[0] + nAbpAddMe[1];
  nMotMme += nAbpAddMe[2];
  // Shift the copied ABPs
  for(n = nAbpCp - 1; n >= 0; n--) {
	UpdateActinAbpMonomerListSubroutine(nAbpMe + nAbpAddSumMe + n,
			nAbpMe + n, 1);
  }
  // Determine the kind of free ABP to add
  for(n = 0; n < nAbpAddSumMe; n++) {
	CS = 0;
	while(CS == 0) {
		kind = GenRandIntIndex(NK_ABP);
		if (nAbpAddMe[kind] > 0) {
			nAbpAddMe[kind]--;
			CS = 1;
		}
	}
	// Add the free ABP
	UpdateActinAbpMonomerListSubroutine2(nAbpMe + n, nAbp + sum + n, kind + 1);
	if ((kind != 2 && gTglImpAcpM == 0) || (kind == 2 && gTglImpMotM == 0)) {
		GenRandPosSubdom(&P2(abp.r,nAbpMe + n,0));
	}
	if (kind != 2 && gTglImpAcpM != 0) {
		InsertElement1dArrayWoChk(acpM.l, &acpM.c, nAbp + sum + n);
	}
	if (kind == 2 && gTglImpMotM != 0) {
		InsertElement1dArrayWoChk(motM.l, &motM.c, nAbp + sum + n);
	}
  }
  // All other free ABPs belonging to other submains should be not here
  for (n = 0; n < nAbpMallSum; n++) {
	if (!(n >= sum && n < sum + nAbpAddSumMe)) {
		iAbp[nAbp + n] = -1;
	}
  }
  nAbpMe += nAbpAddSumMe;
  nAbp += nAbpMallSum;
}

void AddFreeActin(void) {
  int n, k, ind[NDIM], nActAddMe, cnt, *cntAll, sum;
  int begin[NDIM], end[NDIM];
  double dr[NDIM], dist, thres[2], r[NDIM];

  if (nActMall > 0) {
	nActAddMe = nActMall / nCpu;
	sum = rank * nActAddMe;
	for (n = nActCp - 1; n >= 0; n--) {
		UpdateActinAbpMonomerListSubroutine(nActMe + nActAddMe + n,
				nActMe + n, 0);
	}
	for (n = 0; n < nActAddMe; n++) {
		UpdateActinAbpMonomerListSubroutine2(nActMe + n,
				nAct + sum + n, 0);
		actM.l[n + actM.c] = nAct + sum + n;
	}
	for (n = 0; n < nActMall; n++) {
		CONT(n >= sum && n < sum + nActAddMe);
		iAct[nAct + n] = -1;
	}
	actM.c += nActAddMe;
	nActMe += nActAddMe;
	nAct += nActMall;
	if (rank == 0) {
  		Printf0("%d actins are added.\n", nActMall);
	} 
  } 
}

void DeleteAbp(void) {
  int n, k, CS, ind, ind2, side, cnt, nDelSum, *pArr;
  int *abpIndL, *abpKindL, *chAbp2, *abpMotId2, nDel[NK_ABP], nDel2[NK_ABP];
  double *rAbp2;
  ListInt del;

  for(n = 0; n < NK_ABP; n++) {
	if (nAbpMall[n] < 0) {
		nDel[n] = -1 * nAbpMall[n];
		nAbpDet[n] += nAbpMall[n];
		nAbpMall[n] = 0;
	}
	else { nDel[n] = 0; }
  }
  nMot -= nDel[2];

  MALLOC(del.l,int,nAbp);
  MALLOC(chAbp2,int,nAbp*nChAb);
  MALLOC(rAbp2,double,nAbp*NDIM);
  MALLOC(abpIndL,int,nAbp);
  if (motSA.gTgl != 0) {
	MALLOC(abpMotId2,int,nAbp);
  }

  Printf0("%d ACPC, %d ACPB, %d motor are deleted.\n", 
			nDel[0], nDel[1], nDel[2]); 
  nDelSum = V3SUM(nDel);
  del.c = 0;
  if (rank == 0) {
    MALLOC(abpKindL,int,nAbp);
	FOR_ABP(n) { abpKindL[n] = P2A(chAbp,n,2,nChAb); }
	for(n = 0; n < NK_ABP; n++) {
		for(k = 0; k < nDel[n]; k++) {
			CS = 0;
			while (CS == 0) {
				ind = GenRandIntIndex(nAbp);
				CONT(!(abpKindL[ind] == n));
				if (ISMTF(n)) {
					CONT(P2A(chAbp,ind,3,nChAb) > -1 
							&& P2A(chAbp,ind,4,nChAb) > -1);
				}
				del.l[del.c] = ind;
				(del.c)++;
				abpKindL[ind] = -1;
				CS = 1;
			}
		}
	}
	free(abpKindL);
  }
  MPI_Bcast(del.l, nDelSum, MPI_INT, 0, MPI_COMM_WORLD);

  for(n = 0; n < nDelSum; n++) {
	ind = del.l[n];
	for(k = 0; k < 5; k++) {
		CONT(k == 2);
		ind2 = P2A(chAbp,ind,k,nChAb);
		CONT(!(ind2 > -1));
		if (k < 2) {
			side = FindAbpActinChain(ind2, ind, 1);
			P2A(chAct,ind2,side,nChAc) = -1;	
		}
		else {
			P2A(chAbp,ind2,7 - k,nChAb) = -1;	
		}
	}
	SetAllValue1dArrayInt(&P2A(chAbp,ind,0,nChAb), nChAb, -1);
  }

  cnt = 0;
  FOR_ABP(n) {
	CONT(!(P2A(chAbp,n,2,nChAb) > -1));
	Copy1dArrayInt(&P2A(chAbp2,cnt,0,nChAb), &P2A(chAbp,n,0,nChAb), nChAb);
	V3COPY(&P2(rAbp2,cnt,0), &P2(rAbp,n,0));
	if (motSA.gTgl != 0) {
		abpMotId2[cnt] = abpMotId[n];
	}
	for(k = 0; k < 2; k++) {
		ind = P2A(chAbp,n,k,nChAb);
		CONT(!(ind > -1));
		side = FindAbpActinChain(ind, n, 1);
		P2A(chAct,ind,side,nChAc) = cnt;
	}
	abpIndL[n] = cnt;
	cnt++;
  }
  nAbp = cnt;
  FOR_ABP(n) {
	for(k = 3; k < 5; k++) {
		CONT(P2A(chAbp2,n,k,nChAb) < 0);
		P2A(chAbp2,n,k,nChAb) = abpIndL[P2A(chAbp2,n,k,nChAb)];
	}
  }

  FOR_ABP(n) {
    Copy1dArrayInt(&P2A(chAbp,n,0,nChAb), &P2A(chAbp2,n,0,nChAb), nChAb);
    V3COPY(&P2(rAbp,n,0), &P2(rAbp2,n,0));
	CONT(!(motSA.gTgl != 0));
	abpMotId[n] = abpMotId2[n];
  }

  free(chAbp2);
  free(rAbp2);
  free(del.l);
  free(abpIndL);
  if (motSA.gTgl != 0) { free(abpMotId2); }
}

// Elminate free actin segments. Because this code has neither polymerization
// nor depolymerization, they play no role here.
void DeleteFreeActin(void) {
  int n, k, side, cnt, *actInd, abpInd;

  Printf0("\n======================= Eliminating free actin segments "
			"========================\n");
  MALLOC(actInd, int, nAct);

  cnt = 0;
  FOR_ACT(n) {
	CONT(!(P2A(chAct,n,0,nChAc) > -1 || P2A(chAct,n,1,nChAc) > -1));
	actInd[n] = cnt;
	cnt++;
  }
  FOR_ACT(n) {
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(chAct,n,k,nChAc) > -1));
		P2A(chAct,n,k,nChAc) = actInd[P2A(chAct,n,k,nChAc)]; 
	}
	for(k = 2; k < nChAc; k++) {
		CONT(!(P2A(chAct,n,k,nChAc) > -1));
		abpInd = P2A(chAct,n,k,nChAc);
		side = (n == P2A(chAbp,abpInd,0,nChAb)) ? 0 : 1;
		P2A(chAbp,abpInd,side,nChAb) = actInd[n];
	}
  }
  Printf0("%d free actin segments are eliminated.\n", nAct - cnt); 
  if (rheoWay > 0) {
	for(n = 0; n < appStraPar.c; n++) {
		appStraPar.l[n] = actInd[appStraPar.l[n]];
	}
	for(n = 0; n < meaStrePar.c; n++) {
		meaStrePar.l[n] = actInd[meaStrePar.l[n]];
	}
  }
  FOR_ACT(n) {
	CONT(!(P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) < 0));
	for(k = n; k < nAct - 1; k++) { 
		Copy1dArrayInt(&P2A(chAct,k,0,nChAc), &P2A(chAct,k + 1,0,nChAc), 
				nChAc);
		V3COPY(&P2(rAct,k,0), &P2(rAct,k + 1,0));
		actInd[k] = actInd[k + 1];
		fixAct[k] = fixAct[k + 1];
	}
	nAct--;
	n--;
  }
  free(actInd);
  MPI_Barrier(MPI_COMM_WORLD);
}

// Eliminate inactive ABPs if any. This is usally used to delete many inactive
// ABPs in a polymerized network. This should be executed before 
// ExtractConfig(), and it cannot be used at the mid of simulation.
void DeleteInactiveAbp(void) {
  int n, k, cnt, actInd, *chAbp2, side;
  double *rAbp2;

  Printf0("\n========================== Eliminating inactive ABPs "
			"===========================\n");

  MALLOC(chAbp2, int, nAbp*nChAb);
  MALLOC(rAbp2, double, nAbp*3);

  cnt = 0;
  // Information of active ABPs is stored in separate arrays.
  FOR_ABP(n) {
    if (P2A(chAbp,n,0,nChAb) > -1 && P2A(chAbp,n,1,nChAb) > -1) {
		Copy1dArrayInt(&P2A(chAbp2,cnt,0,nChAb), &P2A(chAbp,n,0,nChAb), nChAb);
		V3COPY(&P2(rAbp2,cnt,0), &P2(rAbp,n,0));
        for(k = 0; k < 2; k++) {
			actInd = P2A(chAbp,n,k,nChAb);
			side = FindAbpActinChain(actInd, n, 1);
			P2A(chAct,actInd,side,nChAc) = cnt;
        }
        cnt++;
    }
    else if (P2A(chAbp,n,0,nChAb) > -1 && P2A(chAbp,n,1,nChAb) < 0) {
		actInd = P2A(chAbp,n,0,nChAb);
		side = FindAbpActinChain(actInd, n, 1);
		P2A(chAct,actInd,side,nChAc) = -1;
    }
  }
  Printf0("%d of %d ABPs are deleted since they are in an inactive "
			"state!!\n", nAbp - cnt, nAbp);
  nAbp = cnt;
  nMot = 0;
  // Put it back to original array
  FOR_ABP(n) {
    Copy1dArrayInt(&P2A(chAbp,n,0,nChAb), &P2A(chAbp2,n,0,nChAb), nChAb);
    V3COPY(&P2(rAbp,n,0), &P2(rAbp2,n,0));
	if (P2A(chAbp,n,2,nChAb) == 2) { nMot++; }
  }
  free(chAbp2);
  free(rAbp2);
  MPI_Barrier(MPI_COMM_WORLD);
}

/*----------------------- Adding and deleting elements -----------------------*/

/*------------------- Severing filaments at the beginning --------------------*/

// Sever filaments whose length is greater than the minimum width of the 
// computational domain in order to preclude artifacts due to small domain. 
void SeverLongFilament(void) {
  int n, k, filaLen, curr, next, minDimDom = POS_LARGE_VALUE, ind;

  FOR_NDIM(k) { 
    if ((int)dimDom[k] < minDimDom) { minDimDom = (int)dimDom[k]; }
  }
  FOR_ACT(n) {
	// Find a pointed end
    CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
	filaLen = 2; 
	curr = P2A(chAct,n,0,nChAc);
	while (P2A(chAct,curr,0,nChAc) > -1) {
		curr = P2A(chAct,curr,0,nChAc);
		filaLen++;
		// If filament is longer than minimal width,the filament is severed.
		CONT(!(filaLen >= minDimDom));
		// Go backward by three segments
		for(k = 0; k < 3; k++) { curr = P2A(chAct,curr,1,nChAc); }
		next = P2A(chAct,curr,0,nChAc);
		// Sever it into two filaments
		P2A(chAct,curr,0,nChAc) = -1;
		P2A(chAct,next,1,nChAc) = -1;
		// Detach any ABP on the severed tip
		for(k = 0; k < 2; k++) {    
			ind = (k == 0) ? curr : next;
			DetachAbpOnActin(ind);
		}
		curr = next;
		filaLen = 1;
	}
  }
}

// Subroutine for SeverFilaOnPlane()
void SeverFilaOnPlaneSubroutine(int ind, int direc) {
  fixAct[ind]++;
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC) 
		|| (bulkRheoType == 1 && direc == dirStr))) {
	fixAct[ind] += 10;
	InsertElement1dArrayWChk(meaStrePar.l,&meaStrePar.c,ind);
	InsertElement1dArrayWChk(appStraPar.l,&appStraPar.c,ind);
  }
}

// Sever all filaments crossing boundaries in a designated direction.
// After severing filaments, they are clamped to the boundaries.
// meaStrePar.l and appStraPar.l are filled here.
// mode = 1: clamped, 0: not clamped 
void SeverFilaOnPlane(int direc, int mode) {
  int n, k, l, ind, ind2, cnt, side, CS;
  double longLen;
  char namePlane[3];
 
  FOR_ACT(n) {
	ConvertRectDomainVector(&P2(rAct,n,0), 1);
  }
  FOR_ABP(n) {
	ConvertRectDomainVector(&P2(rAbp,n,0), 1);
  }
  switch(direc) {
    case 0: strcpy(namePlane, "yz"); break;
    case 1: strcpy(namePlane, "zx"); break;
    case 2: strcpy(namePlane, "xy"); break;
  }
  Printf0("\n==================== Severing filaments that cross %s plane "
			"====================\n", namePlane);
  cnt = 0;

  // Delete all ABP bonds near boundaries.
  longLen = 1.0;
  for(n = 0; n < NK_ABP; n++) {
	if (abpF.len[n].n + 0.5 * actF.dia > longLen) {
		longLen = abpF.len[n].n + 0.5 * actF.dia;
	}
  }
  FOR_ABP(n) {
    CONT(!(P2(rAbp,n,direc) >= rGrid[direc][nGrid[direc] - 1] - longLen
	        || P2(rAbp,n,direc) < rGrid[direc][0] + longLen));
	cnt++;
	DetachActinOnAbp(n);
  }
  Printf0("%d ABPs located near %s plane are unbound.\n", cnt, namePlane);

  // Fill in meaStrePar.l and appStraPar.l
  cnt = 0;
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC)
        || (bulkRheoType == 1 && direc == dirStr))) {
	meaStrePar.c = 0;
	appStraPar.c = 0;
  }
  // Find filaments crosslinking boundaries and sever/clamp them.
  FOR_ACT(n) {
	CONT(!(P2(rAct,n,direc) >= rGrid[direc][nGrid[direc] - 1] - 1.0 
			|| P2(rAct,n,direc) <= rGrid[direc][0] + 1.0 ));
	CS = (P2(rAct,n,direc) > rGrid[direc][0] + dimDomH[direc]) ? 1 : 0;
	for(k = 0; k < 2; k++) {
		ind = P2A(chAct,n,k,nChAc);
		CONT(!(ind > -1));
		CONT(!(fabs(P2(rAct,n,direc) - P2(rAct,ind,direc)) > dimDomH[direc]));
		if (P2A(chAct,ind,k,nChAc) > -1
				&& P2A(chAct,n,1 - k,nChAc) > -1) {
			if (mode == 1) {
				for(l = 0; l < 2; l++) {
					ind2 = (l == 0) ? n : ind;
					side = (l == 0) ? 1 - k : k;		
					CONT(!(fixAct[P2A(chAct,ind2,side,nChAc)] < 0));
					fixAct[ind2] = direc * 2 + 1;
					if (rheoWay > 0) {
						fixAct[P2A(chAct,ind2,side,nChAc)] = direc * 2 + 1;
					}
					if ((l == 0 && CS == 1) || (l == 1 && CS == 0)) {
						SeverFilaOnPlaneSubroutine(ind2, direc);
						if (rheoWay > 0) {
							SeverFilaOnPlaneSubroutine(
									P2A(chAct,ind2,side,nChAc), direc);
						}
					}
				}
			}
			if (rheoWay > 0) {
				DetachAbpOnActin(P2A(chAct,n,1 - k, nChAc));
				DetachAbpOnActin(P2A(chAct,ind,k,nChAc));
			}
		}
		else {
			if (rheoWay > 0) { 
				for(l = 0; l < 2; l++) {
					ind2 = (l == 0) ? n : ind;
					DeleteElement1dArray(meaStrePar.l, &meaStrePar.c, ind2);
					DeleteElement1dArray(appStraPar.l, &appStraPar.c, ind2);
				}
			}
		}
		cnt++;
		DetachAbpOnActin(n);
		DetachAbpOnActin(ind);
		P2A(chAct,n,k,nChAc) = -1;
		P2A(chAct,ind,1 - k,nChAc) = -1;
	}
  }
  Printf0("%d filaments crossing %s plane are severed.\n",
			 cnt, namePlane);
  if (rheoWay > 0 && ((bulkRheoType == 0 && direc == dirNoPBC)
        || (bulkRheoType == 1 && direc == dirStr))) {
	qsort(appStraPar.l, appStraPar.c, sizeof(int), CompInt);
	qsort(meaStrePar.l, meaStrePar.c, sizeof(int), CompInt);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/*------------------- Severing filaments at the beginning --------------------*/

/*----------------- Breaking chains between actins and ABPs ------------------*/

// Detach all ABP bonds from the actin segment
void DetachAbpOnActin(int actInd) {
  int n, abpInd;

  for(n = 2; n < nChAc; n++) {
	CONT(!(P2A(chAct,actInd,n,nChAc) > -1));
	abpInd = P2A(chAct,actInd,n,nChAc);
	if (P2A(chAbp,abpInd,0,nChAb) == actInd) {
		P2A(chAbp,abpInd,0,nChAb) = P2A(chAbp,abpInd,1,nChAb);
	}
	P2A(chAbp,abpInd,1,nChAb) = -1;
	P2A(chAct,actInd,n,nChAc) = -1;
  }
}

// Detach all actin bonds from the ABP
void DetachActinOnAbp(int abpInd) {
  int n, actInd, side;

  for(n = 0; n < 2; n++) {
	CONT(!(P2A(chAbp,abpInd,n,nChAb) > -1));
	actInd = P2A(chAbp,abpInd,n,nChAb);
	side = FindAbpActinChain(actInd, abpInd, 1);
	P2A(chAct,actInd,side,nChAc) = -1;
	P2A(chAbp,abpInd,n,nChAb) = -1;
  }
}

/*----------------- Breaking chains between actins and ABPs ------------------*/

/*------------------------ Inspect chain information -------------------------*/

// Check the integrity of loaded network configuration in terms of 
// chain information.
int InspectChainList(void) {
  int n, k, l, CS, CS2, ind1, ind2, ind3;

  CS = 1;
  FOR_ACTME(n) {
	// Check chain information between actin segments in "act.ch"
	ind1 = act.id[n];
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(act.ch,n,k,nChAc) > -1));
		ind2 = P2A(act.ch,n,k,nChAc);
		ind3 = iAct[ind2];
		CONT(!(P2A(act.ch,ind3,1 - k,nChAc) != ind1));
		Printf("Error: chain information between actin segment %d and"
				" %d is wrong at rank %d.\n", ind1, ind2, rank);
		Printf("%d: ", ind1);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,n,l,nChAc));
		}
		Printf("\n%d: ", ind2);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,ind3,l,nChAc));
		}
		Printf("\n\n");
		CS = 0;
	}
	// Check chain information between actin segment and ABP in "act.ch"
	for(k = 2; k < nChAc; k++) {
		CONT(!(P2A(act.ch,n,k,nChAc) > -1));
		ind2 = P2A(act.ch,n,k,nChAc);
		ind3 = iAbp[ind2];
		CONT(!(P2A(abp.ch,ind3,0,nChAb) != ind1 
				&& P2A(abp.ch,ind3,1,nChAb) != ind1));
		Printf("Error: chain information between actin segment %d and"
				" ABP %d is wrong at rank %d.\n", ind1, ind2, rank);
		Printf("%d: ", ind1);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,n,l,nChAc));
		}
		Printf("\n%d: ", ind2);
		for(l = 0; l < nChAb; l++) {
			Printf("%d\t", P2A(abp.ch,ind3,l,nChAb));
		}
		Printf("\n\n");
		P2A(act.ch,n,k,nChAc) = -1;
		CS = 0;
	}
  }
  // Check chain information between actin segment and ABP in "abp.ch"
  FOR_ABPME(n) {
	ind1 = abp.id[n];
	for(k = 0; k < 2; k++) {
		CONT(!(P2A(abp.ch,n,k,nChAb) > -1));
		ind2 = P2A(abp.ch,n,k,nChAb);
		ind3 = iAct[ind2];
		CS2 = FindAbpActinChain(ind3, ind1, 0);
		CONT(!(CS2 == -1));
		Printf("Error: chain information between actin segment %d and"
				" ABP %d is wrong at rank %d.\n", ind2, ind1, rank);
		Printf("%d: ", ind2);
		for(l = 0; l < nChAc; l++) {
			Printf("%d\t", P2A(act.ch,ind3,l,nChAc));
		}
		Printf("\n%d: ", ind1);
		for(l = 0; l < nChAb; l++) {
			Printf("%d\t", P2A(abp.ch,n,l,nChAb));
		}
		Printf("\n\n");
		CS = 0;
	}
  }
  return CS;
}

// Check the integrity of loaded network configuration in terms of chain length
int InspectChainLength(void) {
  int n, k, ind1, ind2, ind3, CS = 1;
  double len, dist;

  // Check chain lengths using "act.ch"
  FOR_ACTME(n) {
	ind1 = act.id[n];	
	for(k = 0; k < nChAc; k++) {
		ind2 = P2A(act.ch,n,k,nChAc);
		CONT(ind2 < 0);
		if (k < 2) { 
			ind3 = iAct[ind2];
			len = CalcDist(&P2(act.r,n,0), &P2(act.r,ind3,0), 0); 
		}
		else { 
			ind3 = iAbp[ind2];
			len = CalcDistActinAbp(ind1, ind2, 0);
		}
		if (k < 2) { dist = AVG2(1., 1.); }
		else { dist = 0.5 * actF.dia + abpF.len[K_ABP(ind3)].n; }

   		CONT(!(len < dist - 0.5 || len > dist + 0.5));
		if (k < 2) {
			Printf("Error: distance between actin segment %d and %d is "
					"wrong at rank %d: %lf (%lf)\n", 
					ind1, ind2, rank, len, dist);
			Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(act.r,n,0), 
					P2(act.r,n,1), P2(act.r,n,2));
			Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(act.r,ind3,0), 
					P2(act.r,ind3,1), P2(act.r,ind3,2));
		}
		else {
			Printf("Error: distance between actin segment %d and ABP"
					" %d is wrong at rank %d: %lf (%lf)\n", 
					ind1, ind2, rank, len, dist);
			Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(act.r,n,0), 
					P2(act.r,n,1),	P2(act.r,n,2));
			Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(abp.r,ind3,0), 
					P2(abp.r,ind3,1), P2(abp.r,ind3,2));
		}
		CS = 0;
	}
  }
  // Check chain lengths using "abp.ch"
  FOR_ABPME(n) {
	ind1 = abp.id[n];
	for(k = 0; k < 2; k ++) {
		ind2 = P2A(abp.ch,n,k,nChAb);
		CONT(ind2 < 0);
		ind3 = iAct[ind2];
		len = CalcDistActinAbp(ind2, ind1, 0);
		dist = 0.5 * actF.dia + abpF.len[K_ABP(n)].n;
		CONT(!(len < dist - 0.5 || len > dist + 0.5));
		Printf("Error: distance between actin segment %d and ABP"
				" %d is wrong at %d: %lf (%lf)\n", ind2, ind1, rank, len, dist);
		Printf("%d: %lf\t%lf\t%lf\n", ind2, P2(act.r,ind3,0), 
				P2(act.r,ind3,1), P2(act.r,ind3,2));
		Printf("%d: %lf\t%lf\t%lf\n", ind1, P2(abp.r,n,0), 
				P2(abp.r,n,1), P2(abp.r,n,2));
		CS = 0;
	}
  }
  return CS;
}

// Inspect chain length and information by calling InspectChainList() and
// InspectChainLength().
void InspectChainListAndLength(void) {
  int CS1, CS2;
  Printf0("\n================ Inspecting the integrity of chain information "
		"=================\n");
  CS1 = InspectChainList();
  CS2 = InspectChainLength();
  if (CS1 != 1) { Printf("Chain list has a problem!!\n"); }
  if (CS2 != 1) { Printf("Chain length has a problem!!\n"); }
  MPI_Barrier(MPI_COMM_WORLD);
}

int BinaryPackActinChainArray(int actInd) {
  int k, ch, cntAbp;

  ch = 0;
  cntAbp = 0;
  for(k = 0; k < nChAc; k++) {
	CONT(!(P2A(act.ch,iAct[actInd],k,nChAc) > -1));
	if (k < 2) { ch += (int)(pow(2, k + 1)); }
	else { cntAbp++; }
  }
  k = 3;
  while(cntAbp != 0) {
	if (cntAbp % 2 == 1) { ch += (int)(pow(2, k)); }
	cntAbp /= 2;
	k++;
  }
  if (gTglActCapAll != 0) { 
	if (act.cap[iAct[actInd]] > -1) {
		ch += 1;
	}
  }
  return ch;
}

/*------------------------ Inspect chain information -------------------------*/

/*--------------------------- Pack information -------------------------------*/

// 0: actin0, 2: stalling at barbed end, 4: stalling by traffic jam; 
// 6: stalling by forces
// 1, 3, 5, 7: same for actin 1
// 8: the kind of ABP
int BinaryPackAbpChainArray(int abpInd) {
  int k, ch, side, actInd, locActInd, locAbpInd, locNextActInd;
  double pMotW, critPMotW;

  critPMotW = 1.0e-10;
  ch = 0;
  locAbpInd = iAbp[abpInd];
  for(k = 0; k < 2; k++) {
	actInd = P2A(abp.ch,locAbpInd,k,nChAb);
	// Check whether any actin is bound on ABP
	if (actInd > -1) { 
		ch += (int)(pow(2, k)); 
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
	}
	else { continue; }
	locNextActInd = iAct[P2A(act.ch,locActInd,0,nChAc)];
	// Check traffic jam
	side = FindAbpActinChain(locActInd, abpInd, 0);
	side = (side - 2) / nChAcY;
	if (side < nChAcX - 1) {
		side = (((side - 2) / nChAcY) + 1) * nChAcY + 2;
		side = FindElementArray(&P2A(act.ch,locActInd,side,nChAc), nChAcY, 
				-1, 0, 1);
	}
	else {
		if (locNextActInd > -1) {
			side = FindElementArray(&P2A(act.ch,locNextActInd,2,nChAc), nChAcY, 
					-1, 0, 1);
			// Check whether ABP is located at the barbed ends
			if (P2A(act.ch,locNextActInd,0,nChAc) < 0) {
				ch += (int)pow(2, k + 2);
			}
		}
	}
	if (side < 0) { 
		ch += (int)pow(2., k + 4);
	}
	CONT(K_ABP(locAbpInd) != 2 || motWalk.gTgl == 0);
	// Check stalling by forces
	pMotW = UpdateMotorWalkingSubroutine(abpInd, k);
	if (pMotW < critPMotW) {
		ch += (int)pow(2., k + 6);
	}
  }
  ch += (int)pow(2., K_ABP(locAbpInd) + 8);
  if (ISMTF(K_ABP(locAbpInd))) {
	if (abp.mId[locAbpInd] > -1) {
		ch += (int)pow(2., abp.mId[locAbpInd] + 11);
	}
  }
  return ch;
}

/*--------------------------- Pack information -------------------------------*/

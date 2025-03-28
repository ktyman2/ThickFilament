// This file contains functions which update things.

/*----------- Update positions of all elements for a next time step ----------*/

// In this function, the positions of actins and ABPs are 
// updated by the explicit Euler equation.
void UpdateNewLocation (void) {
  int m, n, k, ind, ind2, direc, dir, nCh, *pArr;
  double r[NDIM], v[NDIM], disp, dr[3][NDIM], invLen, fac, drag;
  double f[NDIM], mag;

  if (rheoWay > 0) {
	dir = (bulkRheoType == 0) ? dirNoPBC : dirStr;
  }
  // Calculate new positions of actin segments by the Euler equation
  FOR_ACTME(n) {
	CONT(ISACTM(n));
	if (act.fix[n] < 1 || bndReb.gTglSpr != 0) {
		V3SET_ALL(r, 0.);
		FOR_NDIM(k) {
			// Check too large forces of actins.
			if (P2(act.f,n,k) > magUnstF) {
				stopSig = 1;
				RecordError(10);
				RecordErrorTotalForce(act.id[n], 0);
			}
			r[k] = (P2(act.f,n,k) + P2(act.fBr,n,k))* dt;
			// Due to medium velocity without inertia, actin segments
			// are enforced to be displaced.
			if (k == dirStr && rheoWay > 0) {
				r[k] += P2(act.r,n,dir) * stra.dsp / dimDom[dir];
			}
		}
		VV3ADD(&P2(act.r,n,0), r);
		// To avoid sqrt() which is computationally intense.
		disp = V3LEN_SQ(r);
		if (4. * disp > SQR(maxDisp)) {
			maxDisp = 2. * sqrt(disp);
		}
	}
	else if (act.fix[n] > 0 && bnd.gTglActMv != 0 && bndReb.gTglSpr == 0) {
		ind = (act.fix[n] % 10) - 1;
		direc = ind / 2;
		V3SET_ALL(r, 0.);
		FOR_NDIM(k) {
			CONT(k == direc);
			r[k] = (P2(act.f,n,k) + P2(act.fBr,n,k)) 
					/ bnd.drag[ind] * dt;
		}
		VV3ADD(&P2(act.r,n,0), r);
	}
  }
  // Calculate new positions of ABPs by the Euler equation
  FOR_ABPME(n) {
	CONT(ISABPIM(n));
	V3SET_ALL(r, 0.);
	FOR_NDIM(k) {
		// Check too large forces of ABPs
		if (P2(abp.f,n,k) > magUnstF) {
			stopSig = 1;
			RecordError(11);
			RecordErrorTotalForce(abp.id[n], 1);
		}
		r[k] = (P2(abp.f,n,k) + P2(abp.fBr,n,k))
				* abpF.drag[K_ABP(n)].inv * dt;
		if (k == dirStr && rheoWay > 0) {
			r[k] += P2(abp.r,n,dir) * stra.dsp / dimDom[dir];
		}
	}
	VV3ADD(&P2(abp.r,n,0), r);
	// To avoid sqrt() which is computationally intense.
	disp = V3LEN_SQ(r);
	CONT(!(4. * disp > SQR(maxDisp)));
	maxDisp = 2. * sqrt(disp);
  }
}

/*----------- Update positions of all elements for a next time step ----------*/

/*------------------------ Updating a neighboring list -----------------------*/

// Calculate the maximum displacement of particles from reference position,
// and if the maximum is large enough, the neighboring list should be updated.
// This concept is called "Heuristic update".
void MeasureDisplacementForNL(void) {
  double lenSq;
  int n;

  if (maxNeiLenSq != POS_LARGE_VALUE || updSubdSize.tgl == 0) {
	maxNeiLenSq = -1.;
	for(n = 0; n < nActMe; n++) {
		CONT(!(P2A(act.ch,n,0,nChAc) > -1 || P2A(act.ch,n,1,nChAc) > -1));
		lenSq = CalcDist(&P2(act.r,n,0), &P2(act.rPrev,n,0), 1);
	    CONT(!(lenSq > maxNeiLenSq));
		maxNeiLenSq = lenSq; 
	}
  }

  if (maxNeiLenSq > dispHeuSq) {
	UpdateNeighborList();
    UpdatePrevLocationForNL();
  }
  if (maxNeiLenSq == POS_LARGE_VALUE) { maxNeiLenSq = -1; }
}

// Remember the reference points for updating neighboring list.
void UpdatePrevLocationForNL(void) {
  int n;

  for (n = 0; n < nActMe * NDIM; n++) {
    act.rPrev[n] = act.r[n];
  }
}

int UpdateNeighborListSubroutine3(int m, int cInd, int *neighInd, int *kind)  {
  int ind, CS;

  CS = 1;

  *kind = SetKind(cInd, act.cyl.c, nAbpMe + nAbpCp);
  if (*kind == 0) {
	*neighInd = cInd;
  }
  else if (*kind == 1) {
	*neighInd = nAct + abp.id[cInd - act.cyl.c];
  }
  else {
	*neighInd = cInd - act.cyl.c - nAbpMe - nAbpCp;
	if (*neighInd >= memb.unit.c) {
		*neighInd = *neighInd - memb.unit.c + nUnitMb;
	}
 	*neighInd += nAct + nAbp;
  }

  if (m == 0 && *kind == 1) {
	ind = cInd - act.cyl.c;
	if (P2A(abp.ch,ind,1,nChAb) < 0) {
		if ((motReb.gTgl == 0 && K_ABP(ind) == 2)
				|| (acpReb.gTgl == 0 && K_ABP(ind) != 2)) {
			CS = -1;
		}
	}
  }
  return CS;
}
int UpdateNeighborListSubroutine2(int m, int *cInd, int *kind)  {
  int n, k, l, CS2, ind, *pArr[2], dim, idx[2], mode;
  double rPnt[6][NDIM], ratio[2], dr[NDIM], dist, len;

  dim = (dimMbNuc == 3 && (kind[0] == 2 || kind[1] == 2)) ? 3 : 2;
  CS2 = 1;
  if (m == 0) {
	if (kind[0] == 0 && kind[1] == 0) {
		pArr[0] = &P2A(act.cyl.l,cInd[0],0,2);
		pArr[1] = &P2A(act.cyl.l,cInd[1],0,2);
		// If the two actins are adjacent on the same filament
		if ((pArr[0][0] == pArr[1][1]) || (pArr[0][1] == pArr[1][0])) {
			CS2 = -1;
		}
		// If the two actins are connected to the same segment
		if (CS2 == 1) {
			for(k = 0; k < 2; k++) {
				if (P2A(act.ch,iAct[pArr[0][k]],1 - k,nChAc) 
						== pArr[1][1 - k]) {
					CS2 = -1;
					break;
				}
			}
		}
	}
	if (abpF.rep.facStf == 0. && kind[0] + kind[1] == 1 && CS2 == 1) {
		idx[0] = (kind[0] == 0) ? 1 : 0;		
		ind = cInd[idx[0]] - act.cyl.c;
		pArr[0] = &P2A(abp.ch,ind,0,nChAb);
		pArr[1] = &P2A(act.cyl.l,cInd[1 - idx[0]],0,2);
		if (pArr[0][0] == pArr[1][0] || pArr[0][0] == pArr[1][1] 
				|| pArr[0][1] == pArr[1][0] || pArr[0][1] == pArr[1][1]) {
			CS2 = -1;
		}
	}
	// No repulsive forces between ABPs
	if (abpF.rep.facStf == 0. && kind[0] == 1 && kind[1] == 1 && CS2 == 1) {
		CS2 = -1;
	}
  }
  else {
	if (kind[0] < 2 && kind[1] < 2) {
		CS2 = -1;
	}
	// If both of them are membrane segments
	if (kind[0] == 2 && kind[1] == 2 && CS2 == 1) {
		for(k = 0; k < 2; k++) {
			if (cInd[k] < act.cyl.c + nAbpMe + nAbpCp + memb.unit.c) {
				pArr[k] = &P2A(memb.unit.l,cInd[k] - act.cyl.c 
						- nAbpMe - nAbpCp,0,dimMbNuc);
			}
			else {
				pArr[k] = &P2A(memb.unitCp.l,cInd[k] - act.cyl.c 
						- nAbpMe - nAbpCp - memb.unit.c,0,dimMbNuc);
			}
		}
		if (dimMbNuc == 2) {
			if ((pArr[0][0] == pArr[1][1]) || (pArr[0][1] == pArr[1][0])) {
				CS2 = -1;
			}
		}
		else {
			if (memb.idx[iMb[pArr[0][0]]] == memb.idx[iMb[pArr[1][0]]]) {
				CS2 = -1;
			}
		}
	}
  }
  if (CS2 == 1) {
	if ((kind[0] < 2 && kind[1] == 2) || (kind[0] == 1 && kind[1] == 0)) {
		V2SET(idx, 1, 0);
	}
	else {
		V2SET(idx, 0, 1);
	}
	dist = 0.;
	for(n = 0; n < 2; n++) {					
		if (kind[idx[n]] == 0) {
	        for(k = 0; k < 2; k++) {
	            ind = iAct[P2A(act.cyl.l,cInd[idx[n]],k,2)];
	            V3COPY(rPnt[n * dim + k], &P2(act.r,ind,0));
	        }
			dist += actF.dia;
		}
		else if (kind[idx[n]] == 1) {
			ind = cInd[idx[n]] - act.cyl.c;
	        V3COPY(rPnt[n * dim], &P2(abp.r,ind,0));
			dist += abpF.len[K_ABP(ind)].n * 2;
		}
		else {
		    for(k = 0; k < dimMbNuc; k++) {
		        if (cInd[idx[n]] < act.cyl.c + nAbpMe 
						+ nAbpCp + memb.unit.c) {
		            ind = iMb[P2A(memb.unit.l,cInd[idx[n]] - act.cyl.c 
							- nAbpMe - nAbpCp,k,dimMbNuc)];
		        }
	    	    else {
	        	    ind = iMb[P2A(memb.unitCp.l,cInd[idx[n]] - act.cyl.c 
							- nAbpMe - nAbpCp - memb.unit.c,k,dimMbNuc)];
		        }
				V3COPY(rPnt[n * dim + k], &P2(memb.r,ind,0));
			}
			if ((mbSld.abp.gTgl != 0 && kind[idx[1 - n]] == 1) 
					|| (mbSld.act.gTgl != 0 && kind[idx[1 - n]] == 0)) { 
				dist += 2 * mbSld.critDist;
			}
			else {
				dist += memb.thk;
			}
		}
		
    }
	dist *= 0.5;
	if (kind[idx[0]] == 0 && kind[idx[1]] == 0) { mode = 0; }
	else if (kind[idx[0]] == 0 && kind[idx[1]] == 1) { mode = 1; }
	else if (kind[idx[0]] == 1 && kind[idx[1]] == 1) { mode = 2; }
	else if (kind[idx[0]] == 2 && kind[idx[1]] == 1) { mode = 3; }
	else if (kind[idx[0]] == 2 && kind[idx[1]] == 0) { mode = 4; }
	else { mode = 5; }	
    len = CalcRepulsiveForcesSubSubroutine(rPnt, dr, ratio, dist, mode);
    if (len > dist + DT_NL_UPDATE) { CS2 = -1; }
  }
  return CS2;
}

void UpdateNeighborListSubroutine(double r[][NDIM], double *rCen, int mode) {
  int n, k, *pNeiPbc, dim, cnt, CS;

  pNeiPbc = (mode == 0) ? neiPbc : neiPbcMb;
  dim = 2;
  V3SET_ALL(rCen, 0.);
  for(n = 0; n < dim; n++) {
	// First, the sheared domain is converted to a cubical one.
	if (rheoWay > 0 && bulkRheoType == 0) 
	{ ConvertRectDomainVector(r[n], 0); }
	FOR_NDIM(k) {
		// positions of copied particles need to be offset.
		CONT(!(nCell[k] > 1 && pbc[k] == 1));
		if (iCell[k] == 0 && r[n][k] >= rGrid[k][nGrid[k] - 1] - neiEdge) { 
			r[n][k] -= dimDom[k]; 
		}
		else if (iCell[k] == nCell[k] - 1 
				&& r[n][k] < rGrid[k][0] + neiEdge) {
			r[n][k] += dimDom[k]; 
		}
	}
	VV3ADD(rCen, r[n]);
  }
  V3SCALE(rCen, INV(dim));
  if (dim == 3) { 
	FOR_NDIM(k) {
		CONT(!(pNeiPbc[k] == 1));
		CS = 1;
		cnt = 0;
		for(n = 0; n < 3; n++) {
			if (fabs(r[n][k] - r[(n + 1) % 3][k]) > dimDomH[k]) { CS = 0; }
			if (r[n][k] > dimDomH[k]) { cnt++; }
		}
		CONT(CS == 1);
		if (rCen[k] > (double)cnt * dimDom[k] / 3.) { 
			rCen[k] -= (double)cnt * dimDom[k] / 3;
		}
		else {
			rCen[k] += (double)(3 - cnt) * dimDom[k] / 3;
		}
	}
  }
  else {
	// If the number of cell in a direction is one with PBC, 
	// their positions need to be adjusted.
	FOR_NDIM(k) {
		CONT(!(pNeiPbc[k] == 1));
		CONT(!(fabs(r[0][k] - r[1][k]) > dimDomH[k]));
		rCen[k] += dimDomH[k] * ((rCen[k] >= rGrid[k][0] + dimDomH[k]) 
				? -1. : 1.);
	}
  }
}

// Update neighboring list for calculating repulsive forces. Due to this list,
// it is not needed to check all distances between particles.
void UpdateNeighborList(void) {
  int m, n, k, l, CS, CS2, offset, nCh, cel, mbCylC;
  int cInd[2], cellPos[2][NDIM], cellInd[2], neighInd[2];
  int ind[3], idCell[NDIM], *pI, *pNeiPbc, kind[2];
  int oft[14][3] = {0,0,0,1,0,0,1,1,0,0,1,0,-1,1,0,0,0,1,1,0,1,1,1,1,0,1,1,
		-1,1,1,-1,0,1,-1,-1,1,0,-1,1,1,-1,1};
  int tOffset, vOffList[][14] = OFFSET_LIST, vOffTableLen[] = OFFSET_LEN, iof;
  double rPnt[3][NDIM], rCen[NDIM]; 
  Cell *pL;
  ListInt *pL2;

  MALLOC(cell.l,int,nActMin+nAbpMin+V3PROD(cell.n));
  memset(cell.l, -1, (nActMin + nAbpMin + V3PROD(cell.n)) * sizeof(int));
  // Build the list of cylindrical segments for actins
  act.cyl.c = 0;
  FOR_ACTMECP(n) {
	ind[0] = P2A(act.ch,n,0,nChAc);
	CONT(!(ind[0] > -1));
	CONT(!(iAct[ind[0]] > -1 && iAct[ind[0]] < nActMe + nActCp));
	V2SET(&P2A(act.cyl.l,act.cyl.c,0,2), act.id[n], ind[0]);
	(act.cyl.c)++;
  }
  memb.unitCp.c = 0;
  if (dimMbNuc == 2) {
	FOR_MBCP(n) {
		ind[0] = P2A(memb.ch,nMbMe + n,0,nChMb);
		CONT(!(iMb[ind[0]] >= nMbMe && iMb[ind[0]] < nMbMe + nMbCp));
		V2SET(&P2A(memb.unitCp.l,memb.unitCp.c,0,2), 
				memb.id[nMbMe + n], ind[0]);
		(memb.unitCp.c)++;
	}
  }
  else {
	FOR_MBCP(n) {
		pI = &P2A(memb.ch,nMbMe + n,0,nChMb);
		nCh = (pI[5] < 0) ? 5 : 6;
		for(k = 0; k < nCh; k++) {
			for(l = 0; l < 2; l++) {
				BREAK(pI[(k + l) % nCh] < memb.id[nMbMe + n]);
				BREAK(!(iMb[pI[(k + l) % nCh]] >= nMbMe 
					&& iMb[pI[(k + l) % nCh]] < nMbMe + nMbCp));
			}
			CONT(l != 2);
			V3SET(&P2A(memb.unitCp.l,memb.unitCp.c,0,3), memb.id[nMbMe + n],
					pI[(k + 1) % nCh], pI[k]);
			(memb.unitCp.c)++;
		}
	}
  }
  if (nAct + nAbp > 0) {
	CheckArraySize(&neigh, &neigh.siz, 2, 0);
	neigh.c = 0;
  }

  for(n = 0; n < act.cyl.c + nAbpMe + nAbpCp + memb.unit.c 
		+ memb.unitCp.c; n++) {

	kind[0] = SetKind(n, act.cyl.c, nAbpMe + nAbpCp);
    if (kind[0] == 0) {
		for(k = 0; k < 2; k++) {
			ind[k] = iAct[P2A(act.cyl.l,n,k,2)];
	        V3COPY(rPnt[k], &P2(act.r,ind[k],0));
		}
    }
    else if (kind[0] == 1) { 
		CONT(abpF.rep.facStf == 0. && n >= act.cyl.c + nAbpMe);
		ind[0] = n - act.cyl.c;
		pI = &P2A(abp.ch,ind[0],0,nChAb);
		CONT(ISABPIM(ind[0]));
       	V3COPY(rPnt[0], &P2(abp.r,ind[0],0));
       	V3COPY(rPnt[1], &P2(abp.r,ind[0],0));
    }
	else {
		for(k = 0; k < dimMbNuc; k++) {
			if (n < act.cyl.c + nAbpMe + nAbpCp + memb.unit.c) {
				ind[k] = iMb[P2A(memb.unit.l,n - act.cyl.c 
						- nAbpMe - nAbpCp,k,dimMbNuc)];
			}
			else {
				ind[k] = iMb[P2A(memb.unitCp.l,n - act.cyl.c - nAbpMe 
						- nAbpCp - memb.unit.c,k,dimMbNuc)];
			}
			BREAK(ind[k] < 0 || ind[k] >= nMbMe + nMbCp);
	        V3COPY(rPnt[k], &P2(memb.r,ind[k],0));
		}
		CONT(k < dimMbNuc);
	}

	for(k = 0; k < 2; k++) {
		if (k == 0) {
			if (kind[0] == 0) { 
				CONT(actF.rep.facStf == 0. && abpF.rep.facStf == 0.
						&& motReb.gTgl == 0 && acpReb.gTgl == 0);
			}
			else if (kind[0] == 1) {
				if (abpF.rep.facStf == 0.) {
					if (pI[1] < 0) {
						CONT((motReb.gTgl == 0 && K_ABP(ind[0]) == 2) 
								|| (acpReb.gTgl == 0 && K_ABP(ind[0]) != 2));
					}
					else { continue; }
				}
			}
			else { continue; }
		}
		else { continue; }
		pL = (k == 0) ? &cell : &cellMb;
		UpdateNeighborListSubroutine(rPnt, rCen, 
				((k == 0) ? 0 : ((kind[0] == 2) ? 2 : 1)));
		if (k == 1 && dimMbNuc == 2) {
			rCen[dirNormMbNuc] = dimDomH[dirNormMbNuc];
		}
		VV3SUB(rCen, pL->base);
		VV3DIV(rCen, pL->wid);
		FOR_NDIM(l) {
			idCell[l] = (int)rCen[l];
			if (idCell[l] >= pL->n[l]) { idCell[l] = pL->n[l] - 1; }
			else if (idCell[l] < 0) { idCell[l] = 0; }
		}
		V3IND_BACK_INT(cel, idCell, pL->n);
		cel += act.cyl.c + nAbpMe + nAbpCp
				+ ((k == 0) ? 0 : memb.unit.c + memb.unitCp.c);
		pL->l[n] = pL->l[cel];
		pL->l[cel] = n;
	}
  }  
  for(m = 0; m < 2; m++) {
	if (m == 0) {
		CONT(nAct + nAbp == 0);
		pL = &cell;
		pL2 = &neigh;
		pNeiPbc = neiPbc;
		mbCylC = 0;
	}
	else {
		continue;
		pL = &cellMb;
		pL2 = &neighMb;
		pNeiPbc = neiPbcMb;
		mbCylC = memb.unit.c + memb.unitCp.c;
	}
	for (n = 0; n < V3PROD(pL->n); n++) {
		V3IND_ASSIGN_INT(n, pL->n, 1, cellPos[0]);
	    cellInd[0] = n + act.cyl.c + nAbpMe + nAbpCp + mbCylC;

		tOffset = 13;
	    if (cellPos[0][2] == pL->n[2] - 1 && pNeiPbc[2] == 0) 
		{ tOffset -= 9; }
		if (pNeiPbc[1] == 0) {
			if (cellPos[0][1] == 0) { tOffset -= 3; }
			else if (cellPos[0][1] == pL->n[1] - 1) { tOffset += 3; }
		}
		if (pNeiPbc[0] == 0) {
			if (cellPos[0][0] == 0) { tOffset -= 1; }
			else if (cellPos[0][0] == pL->n[0] - 1) { tOffset += 1; }
		}
	    cInd[0] = pL->l[cellInd[0]];
	    while (cInd[0] > -1) {
			CS = UpdateNeighborListSubroutine3(m, cInd[0], 
					&neighInd[0], &kind[0]);
			if (CS != 1) {
	        	cInd[0] = pL->l[cInd[0]];
				continue;
			}
	        for (offset = 0; offset < vOffTableLen[tOffset]; offset++) {
	            CS = 1; 
				cInd[1] = -1;
				V3ADD(cellPos[1], cellPos[0], oft[vOffList[tOffset][offset]]);
				FOR_NDIM(k) {
					CONT(!(pNeiPbc[k] == 1));
					if (cellPos[1][k] < 0) { 
						cellPos[1][k] = pL->n[k] - 1; 
					}
					else if (cellPos[1][k] >= pL->n[k]) { 
						cellPos[1][k] = 0;
					}
				}
	            if (CS == 1) {
					cellInd[1] = V3IND_BACK_INT(cellInd[1], cellPos[1], pL->n)
							+ act.cyl.c + nAbpMe + nAbpCp + mbCylC;
	                cInd[1] = pL->l[cellInd[1]];
	            }
	            while (cInd[1] > -1 && CS == 1) {
					CS = UpdateNeighborListSubroutine3(m, cInd[1], 
								&neighInd[1], &kind[1]);
					if (CS != 1) {
			        	cInd[1] = pL->l[cInd[1]];
						continue;	
					}
	                CS2 = 1;
	 				if (!(cellInd[0] != cellInd[1] || 
							(cellInd[0] == cellInd[1] && cInd[1] < cInd[0]))) 
					{ CS2 = -1; }
					if (CS2 == 1) {
						CS2 = UpdateNeighborListSubroutine2(m, cInd, kind);
						if (CS2 == 1) {
							if ((kind[0] == 1 && kind[1] % 2 == 0) 
									|| (kind[0] < 2 && kind[1] == 2)) {
								V2SET(&P2A(pL2->l,pL2->c,0,2), 
										neighInd[1], neighInd[0]);
							}
							else {
								V2SET(&P2A(pL2->l,pL2->c,0,2), 
										neighInd[0], neighInd[1]);
							}
							(pL2->c)++;
						}
	                }
	                cInd[1] = pL->l[cInd[1]];
	            } 
	        } 
	        cInd[0] = pL->l[cInd[0]];
		} 
	}
  }
  free(cell.l);
}

void DeleteActinSegmentInNeighborList(int *actInd) {
  int m, n, ind;  
  ListInt *pL;

  ind = Find2ElementArray(act.cyl.l, act.cyl.c, actInd[0], actInd[1], 0, 2);
  if (ind > -1) {
	DeleteElementArrayByIndex(act.cyl.l, &act.cyl.c, ind, 2);
	for(m = 0; m < 2; m++) { 
		CONT(m == 1);
		pL = (m == 0) ? &neigh : &neighMb;
		for(n = 0; n < pL->c * 2; n++) {
			if (pL->l[n] == ind) {
				DeleteElementArrayByIndex(pL->l, &pL->c,(int)(n / 2), 2);
				n -= n % 2 + 1;
			}
			else if (pL->l[n] > ind && pL->l[n] < nAct) {
				(pL->l[n])--;
			}
		}
	}
  }
}

// kind = 0: actin, 1: ABP, 2: membrane
void DeleteElementInNeighborList(int ind, int kind) {
  int m, n, neighInd, max;
  ListInt *pL;

  if (kind == 0) {
	neighInd = FindElementArray(act.cyl.l, act.cyl.c * 2, ind, 0, 1);
	if (neighInd > -1) {
		neighInd = (int)(neighInd / 2);
		DeleteElementArrayByIndex(act.cyl.l, &act.cyl.c, neighInd, 2);
	}
	max = nAct;
  }
  else if (kind == 1) {
	neighInd = ind + nAct;
	max = -1;
  }
  else if (kind == 2) {
	neighInd = ind + nAct + nAbp; 
	max = nAct + nAbp + nUnitMb;
  }
  else {
	neighInd = ind + nAct + nAbp + nUnitMb; 
	max = nAct + nAbp + nUnitMb * 2; 
  }
  if (neighInd > -1) {
	for(m = 0; m < 2; m++) { 
		CONT(m == 0 && kind == 2);
		CONT(m == 1);
		pL = (m == 0) ? &neigh : &neighMb;
		for(n = 0; n < pL->c * 2; n++) {
			if (pL->l[n] == neighInd) {	
				DeleteElementArrayByIndex(pL->l, &pL->c,(int)(n / 2),2);
				n -= n % 2 + 1;
			}
			else if (pL->l[n] > neighInd && pL->l[n] < max) {	
				(pL->l[n])--;
			}
		}
	}
  }
}

// mode = 0~1: actin, 2: ABP
void InsertElementInNeighborList(int ind1, int ind2, int mode) {
  int n, k, CS, neighInd, locInd, idx, ind3, sft, mode2, kind;
  int oppActInd[2], *cylInd, *pArr;
  double len, dist, dr[NDIM], ratio[4], rPnt[3][NDIM], rPnt2[6][NDIM];
  ListInt *pL;
  // actin
  if (mode < 2) {
	for(n = 0; n < 2; n++) {
		V3COPY(rPnt[n],&P2(act.r,iAct[(n == 0) ? ind2 : ind1],0));
		oppActInd[n] = P2A(act.ch,iAct[(n == 0) ? ind1 : ind2],
				1 - (mode + n) % 2,nChAc);
	}
	neighInd = act.cyl.c;
	if (actF.rep.facStf > 0.) {
		for(n = 0; n < 2; n++) {
			V3COPY(rPnt2[n], rPnt[n]);
		}
		sft = 2;
		mode2 = 0;
		dist = actF.dia;
		CS = 1;
	}
	else { CS = 0; }
  }
  // ABP
  else if (mode == 2) {
	locInd = iAbp[ind1];
	V3COPY(rPnt[0],&P2(abp.r,locInd,0));
	neighInd = nAct + ind1;	
	kind = K_ABP(locInd);
	CS = 1;
	if (P2A(abp.ch,locInd,1,nChAb) < 0) {
		if (abpF.rep.facStf == 0. && ((acpReb.gTgl == 0 && kind != 2) 
				|| (motReb.gTgl == 0 && kind == 2))) {
			CS = 0;
		}
	}
	else {
		if (abpF.rep.facStf == 0.) {
			CS = 0;
		}
	}
	if (CS == 1) {
		V3COPY(rPnt2[2], rPnt[0]);
		sft = 0;
		mode2 = 1;
		dist = AVG2(actF.dia, abpF.len[kind].n * 2);
	}
  }
  // membrane
  else {
	pArr = (mode == 3) ? &P2A(memb.unit.l,ind1,0,dimMbNuc) 
			: &P2A(memb.unitCp.l,ind1,0,dimMbNuc);
	if (iMb[pArr[0]] < 0 || iMb[pArr[1]] < 0) { return; }
	if (dimMbNuc == 3) { 
		if (iMb[pArr[2]] < 0) { return; } 
	}
	for(n = 0; n < dimMbNuc; n++) {
		V3COPY(rPnt[n],&P2(memb.r,iMb[pArr[n]],0));
	}	
	idx = memb.idx[iMb[pArr[0]]];
	neighInd = nAct + nAbp + ind1 + ((mode == 3) ? 0 : nUnitMb);

	for(n = 0; n < dimMbNuc; n++) {
		V3COPY(rPnt2[n], rPnt[n]);
	}
	sft = dimMbNuc;
	mode2 = 4;
	dist = AVG2(actF.dia, memb.thk);
	CS = 1;
  }

  if (CS == 1) {
	for(n = 0; n < act.cyl.c; n++) {
		cylInd = &P2A(act.cyl.l,n,0,2);
		// If two actin segments are connected, avoid this pair.
		if (mode < 2) {
			for(k = 0; k < 2; k++) {
				ind3 = (k == 0) ? ind1 : ind2;
				CONT(cylInd[0] == ind3 || cylInd[0] == oppActInd[k]
						|| cylInd[1] == ind3 || cylInd[1] == oppActInd[k]);
			}
		}
		CONT(iAct[cylInd[0]] < 0 || iAct[cylInd[1]] < 0);
		V3COPY(rPnt2[sft],&P2(act.r,iAct[cylInd[0]],0));
		V3COPY(rPnt2[1 + sft],&P2(act.r,iAct[cylInd[1]],0));

	    len = CalcRepulsiveForcesSubSubroutine(rPnt2, dr, ratio, dist, mode2);
		CONT(!(len < dist + DT_NL_UPDATE));
		if (mode >= 3) {
			V2SET(&P2A(neighMb.l,neighMb.c,0,2), neighInd, n);
			(neighMb.c)++;
		}
		else {
			V2SET(&P2A(neigh.l,neigh.c,0,2), n, neighInd);
			(neigh.c)++;
		}
	}
  }

  if (mode < 2) {
	CS = (abpF.rep.facStf > 0. || acpReb.gTgl != 0 || motReb.gTgl != 0) ? 1 : 0;
	if (CS == 1) {
		for(n = 0; n < 2; n++) {
			V3COPY(rPnt2[n], rPnt[n]);
		}
		sft = 2;
		mode2 = 1;
		dist = 0.5 * actF.dia;
	}
  }
  else if (mode == 2) {
	CS = (abpF.rep.facStf > 0.) ? 1 : 0;
	if (CS == 1) {
		V3COPY(rPnt2[2], rPnt[0]);
		sft = 0;
		dist = 0.5 * abpF.len[K_ABP(locInd)].n;
	}
	mode2 = 2;
  }
  else {
	CS = 1;
	for(n = 0; n < dimMbNuc; n++) {
		V3COPY(rPnt2[n], rPnt[n]);
	}
	sft = dimMbNuc;
	mode2 = 3;
	dist = (mbSld.gTgl != 0) ? mbSld.critDist : 0.5 * memb.thk;
  }
  if (CS == 1) {
	FOR_ABPMECP(n) {
		pArr = &P2A(abp.ch,n,0,nChAb);
		kind = pArr[2];
		CONT(abp.id[n] == ind1);
		CONT(ISABPIM(n));
		if (mode < 2) {	
			if (pArr[1] < 0) {
				CONT(abpF.rep.facStf == 0. && ((acpReb.gTgl == 0 && kind != 2) 
						|| (motReb.gTgl == 0 && kind == 2)));
			}
			else {
				CONT(abpF.rep.facStf == 0.); 
			}
		}
		else if (mode == 2) {
			if (motSA.gTgl != 0 && K_ABP(locInd) == 2 && kind == 2) {
				CONT(pArr[3] == ind1 || pArr[4] == ind1);
			}
		}
		V3COPY(rPnt2[sft], &P2(abp.r,n,0));
	   	len = CalcRepulsiveForcesSubSubroutine(rPnt2, dr, ratio, 
				dist, mode2);
		CONT(!(len < dist + DT_NL_UPDATE + 0.5 * abpF.len[kind].n));
		if (mode < 2) {
			V2SET(&P2A(neigh.l,neigh.c,0,2), act.cyl.c, nAct + abp.id[n]); 
			(neigh.c)++;
	 	}
		else if (mode == 2) {
			V2SET(&P2A(neigh.l,neigh.c,0,2), neighInd, nAct + abp.id[n]); 
			(neigh.c)++;
		}
		else {
			V2SET(&P2A(neighMb.l,neighMb.c,0,2), neighInd, nAct + abp.id[n]); 
			(neighMb.c)++;
		}
	}
  }

  if (mode < 2) {
	P2A(act.cyl.l,act.cyl.c,mode,2) = ind1;
	P2A(act.cyl.l,act.cyl.c,1 - mode,2) = ind2;
	(act.cyl.c)++;
  }
}

/*------------------------ Updating a neighboring list -----------------------*/

/*------------------------ Dynamic behaviors of actins -----------------------*/

void UpdateActinSeveringSubroutine(int *actInd) {
}
 
void UpdateActinSeverAnnealSubroutine(int actInd, int iFila,
		ListInt *sendL, int mode) {
}

void UpdateActinSevering(void) {
}

void UpdateActinNucleation(void) {
  int n, k, CS, cnt, nRep, actInd[2], locActInd[2], ind[2];
  double pActNuc, dr[NDIM], fac, dimDomC[NDIM];

  V3SUB(dimDomC, &P2A(bnd.r,1,0,NDIM), &P2A(bnd.r,0,0,NDIM));
  if (actNuc.gTglFN == 0) {
	fac = (rheoWay > 0 && bulkRheoType == 1) ? (dimDom[dirStr] 
			/ (P2A(rGridInit,1,dirStr,NDIM) - P2A(rGridInit,0,dirStr,NDIM))) 
			: 1.;
	fac /= V3PROD(dimDomC);
	pActNuc = 1. - exp(actNuc.facP * (double)actM.c * fac);
	nRep = actM.c * nActPerSeg;
  }
  else { 
	pActNuc = 1.; 
	nRep = actNuc.cntFNme;
  }
  pActNuc = AdjustDynamicsRate(pActNuc);
  for(n = 0; n < nRep; n++) {
	BREAK(actM.c < 2);
	CONT(!(genrand_real3() < pActNuc));
	cnt = 0;
	// Choose actin monomers from the list
	for(k = 0; k < actM.c; k++) {
		CS = FindElementArray(noActDyn.l, noActDyn.c, actM.l[k], 0, 2);
		CONT(CS > -1);
		ind[cnt] = k;
		cnt++;
		BREAK(cnt == 2);
	}
    // If there are not two available actin monomers, there is no reason to
    // continue the for loop.
	BREAK(!(cnt == 2));
	// Delte actin monomers from the list
	for(k = 1; k >= 0; k--) {
		actInd[k] = actM.l[ind[k]];
		DeleteElementArrayByIndex(actM.l, &actM.c, ind[k], 1);
	}
	pActNuc = 1. - exp(actNuc.facP * (double)actM.c * fac);
	for(k = 0; k < 2; k++) {
		locActInd[k] = iAct[actInd[k]];
		act.len[locActInd[k]] = k;
		P2A(act.ch,locActInd[k],k,nChAc) = actInd[1 - k];
		act.iF[locActInd[k]] = iFilaP.l[0];
		act.fix[locActInd[k]] = -1;
		if (bndMat.gTgl != 0 && bndUnb.gTgl != 0)
		{ act.nFA[locActInd[k]] = 0; }
		if (gTglActCapAll != 0) { act.cap[locActInd[k]] = -1; }
		if (actAge.gTgl != 0) { act.age[locActInd[k]] = currTimeStep; }
		V3SET_ALL(&P2(act.fBr,locActInd[k],0), 0.);
	}
	DeleteElementArrayByIndex(iFilaP.l, &iFilaP.c, 0, 1);
	CS = 0;
	// If it is out of boundary, nucleation cannot happen
	while (CS != 1) {
		CS = 1;
		FOR_NDIM(k) {
			P2(act.r,locActInd[0],k) = P2A(bnd.r,0,k,NDIM) 
					+ dimDomC[k] * genrand_real3();
		}
		GenRandDirecVec(dr);
		V3ADD(&P2(act.r,locActInd[1],0), &P2(act.r,locActInd[0],0), dr);
		ApplyBoundCondVector(&P2(act.r,locActInd[1],0), -1, 0);
		CS = CheckParticleInDomain(&P2(act.r,locActInd[1],0));
		CONT(CS != 1);
		if (bnd.gTglRnd != 0) {
			CS = CheckActinAbpOverlapBoundary(&P2(act.r,locActInd[1],0));
			CONT(CS != 1);
		}
	}
	for(k = 0; k < 2; k++) {
		V3COPY(&P2(act.rPrev,locActInd[k],0), &P2(act.r,locActInd[k],0));
	}
	InsertElementInNeighborList(actInd[0], actInd[1], 0);
	(actNuc.cntMe)++;
	nActFilaMe++;
  	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), actInd[0], currTimeStep);
  	V2SET(&P2A(noActDyn.l,noActDyn.c + 1,0,2), actInd[1], currTimeStep);
  	(noActDyn.c) += 2;
	if (actNuc.gTglFN != 0) { actNuc.cntFNme--;	}
	cntNucAssAnn++;
  }
}

void UpdateActinAssembly(void) {
  int n, k, side, CS, actInd, locActInd, *pArr, chkBst;
  double pActAss[2], r[NDIM], dr[NDIM], dimDomC[NDIM], fac;

  fac = (rheoWay > 0 && bulkRheoType == 1) ? (dimDom[dirStr] 
		/ (P2A(rGridInit,1,dirStr,NDIM) - P2A(rGridInit,0,dirStr,NDIM))) : 1.;
  V3SUB(dimDomC, &P2A(bnd.r,1,0,NDIM), &P2A(bnd.r,0,0,NDIM));
  fac /= V3PROD(dimDomC);
  for(n = 0; n < 2; n++) {
	pActAss[n] = 1. - exp(actAss.facP[n] * (double)actM.c * fac);
	pActAss[n] = AdjustDynamicsRate(pActAss[n]);
  }
  FOR_ACTME(n) {
	CONT(ISACTM(n));

	BREAK(actM.c == 0);
	// If an actin is subject to the bursting depolymerization, 
	// polymerization shouldn't occur on it.
	chkBst = -1;
	if (actBst.tgl != 0) {
		chkBst = FindElementArray(actBst.fil.l, actBst.fil.c, 
				act.iF[n], 0, 3);
	}
	CONT(chkBst > -1);

	pArr = &P2A(act.ch,n,0,nChAc);
	// If barbed or pointed end
	CONT(!((pArr[1] > -1 && pArr[0] < 0) || (pArr[0] > -1 && pArr[1] < 0)));
	if (gTglActCapAll != 0) {
		CONT(act.cap[n] > -1);
	}
	// Check whether a dynamic behavior is allowed on this actin
	CS = FindElementArray(noActDyn.l, noActDyn.c, act.id[n], 0, 2);
	CONT(CS > -1);
	side = (pArr[1] > -1 && pArr[0] < 0) ? 0 : 1;
	// Calculate the position of the particle
	CalcUnitVec(dr,&P2(act.r,n,0), &P2(act.r,iAct[pArr[1 - side]],0));
	V3ADD(r,&P2(act.r,n,0),dr);
	if (bnd.gTglRnd != 0) {
		CS = CheckActinAbpOverlapBoundary(r);
		CONT(CS != 1);
	}
	ApplyBoundCondVector(r, -1, 0);
	// Check probability of assembly
	CONT(!(genrand_real3() < pActAss[side]));
    // Update chain information
    for(k = 0; k < actM.c; k++) {
        CS = FindElementArray(noActDyn.l, noActDyn.c, actM.l[k], 0, 2);
        CONT(CS > -1);
        break;
    }
    CONT(k == actM.c && CS > -1);
    actInd = actM.l[k];
    DeleteElementArrayByIndex(actM.l, &actM.c, k, 1);
	for(k = 0; k < 2; k++) {
		pActAss[k] = 1. - exp(actAss.facP[k] * (double)actM.c * fac);
		pActAss[k] = AdjustDynamicsRate(pActAss[k]);
	}

	locActInd = iAct[actInd];

	act.len[locActInd] = act.len[n] + 1;

	V3COPY(&P2(act.r,locActInd,0),r);
	V3SET_ALL(&P2(act.fBr,locActInd,0), 0.);
	pArr[side] = actInd;
	P2A(act.ch,locActInd,1 - side,nChAc) = act.id[n];
	act.fix[locActInd] = -1;
	if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) { act.nFA[locActInd] = 0; }
	if (gTglActCapAll != 0) { act.cap[locActInd] = -1; }
	if (actAge.gTgl != 0) { act.age[locActInd] = currTimeStep; }
	act.iF[locActInd] = act.iF[n];
	(actAss.cntMe)++;
	V3COPY(&P2(act.rPrev,locActInd,0), &P2(act.r,locActInd,0));
	InsertElementInNeighborList(act.id[n], actInd, side);
	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), act.id[n], currTimeStep);
	V2SET(&P2A(noActDyn.l,noActDyn.c + 1,0,2), actInd, currTimeStep);
	(noActDyn.c) += 2;
	cntNucAssAnn++;
  }
}

void UpdateActinDisassemblySubroutine(int actInd, int side) {
  int locActInd, mbInd, sideMb;

  locActInd = iAct[actInd];
  if (locActInd > -1) {
	P2A(act.ch,locActInd,side,nChAc) = -1;
	// If the remaining part is a monomer, the number of filaments is decreased
	// by one, and additional procedures are performed.
	if (P2A(act.ch,locActInd,1 - side,nChAc) < 0) {
		V3SET_ALL(&P2(act.r,locActInd,0), 0.);
		SetAllValue1dArrayInt(&P2A(act.ch,locActInd,0,nChAc), nChAc, -1);
		if (locActInd < nActMe) {
			act.fix[locActInd] = -1;
			if (bndMat.gTgl != 0 && bndUnb.gTgl != 0)
			{ act.nFA[locActInd] = 0; }
			if (gTglActCapAll != 0) { act.cap[locActInd] = -1; }
			if (actAge.gTgl != 0) { act.age[locActInd] = 0; }
			if (actSev.gTgl != 0 || actNuc.gTgl != 0) {
				InsertElementArrayByIndex(iFilaP.l, &iFilaP.c, 
						&act.iF[locActInd], 0, 1);
			}
			InsertElement1dArrayWoChk(actM.l, &actM.c, actInd);
			if (recTraj.tgl != 0) {
				ReplaceElementInTrajList(actInd, 0);
			}
			(actDis.cntMe)++;
			nActFilaMe--;
		}
		act.iF[locActInd] = -1;
  		if (actNuc.gTglFN != 0 && locActInd < nActMe) { actNuc.cntFNme++; }
	}
  	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), actInd, currTimeStep);
  	(noActDyn.c)++;
  }
}

void UpdateActinDisassemblySubroutine2(int abpInd, int actInd) {
  int abpSide, locAbpInd, *pArr;

  locAbpInd = iAbp[abpInd];
  if (locAbpInd > -1) {
	pArr = &P2A(abp.ch,locAbpInd,0,nChAb);
	abpSide = (pArr[0] == actInd) ? 0 : 1;
	// If active ABP
	if (pArr[1 - abpSide] > -1) {
		UpdateActiveAbpUnbindingSubroutine(actInd, abpInd);
		if (locAbpInd < nAbpMe) {
			if (tglNeiAbpSC != 0 && tglNeiAbpDC == 0) {
				InsertElementInNeighborList(abpInd, -1, 2);
			}
		}
	}
	// If inactive ABP
	else {
		if (tglNeiAbpSC != 0) {
			if (locAbpInd < nAbpMe) {
				if ((K_ABP(locAbpInd) == 2 && gTglImpMotM != 0 
						&& motSA.gTgl == 0) || (K_ABP(locAbpInd) != 2 
						&& gTglImpAcpM != 0) || (K_ABP(locAbpInd) == 2 
						&& motSA.gTgl != 0 && pArr[3] < 0 && pArr[4] < 0)) {
					DeleteElementInNeighborList(abpInd, 1);
				}
			}
		}
		UpdateInactAbpUnbindingSubroutine(actInd, abpInd);
	}
  }
}

void UpdateActinDisassembly(void) {
  int n, k, side, sideMb, *pArr, *pArr2, CS, chkBst, cnt, end;
  int ind, actInd, locActInd, actInd2, locActInd2, abpInd;
  int mbInd, locMbInd;
  double pDis;

  end = (actDgd.gTgl != 0) ? actDgd.c : nActMe;
  for(n = 0; n < end; n++) {
	actInd = (actDgd.gTgl != 0) ? actDgd.l[n] : act.id[n];
	locActInd = iAct[actInd];
	CONT(ISACTM(locActInd));
	CONT(!(locActInd < nActMe));
	pArr = &P2A(act.ch,locActInd,0,nChAc);
	// If it is not one of the two ends, or if the actin is fixed 
	// for some reason, it cannot be disassembled
	CONT(pArr[0] > -1 && pArr[1] > -1);
	side = (pArr[1] > -1 && pArr[0] < 0) ? 0 : 1;
	actInd2 = pArr[1 - side];
	locActInd2 = iAct[actInd2];
	pArr2 = &P2A(act.ch,locActInd2,0,nChAc);

	CS = -1;
	for(k = 0; k < noActDyn.c; k++) {
  		CONT(!(P2A(noActDyn.l,k,0,2) == actInd 
				|| P2A(noActDyn.l,k,0,2) == actInd2));
		CS = 1;
		break;
	}
	CONT(CS != -1);
	cnt = HowManyAbpActinChain((side == 0) ? locActInd2 : locActInd, 0);
	chkBst = -1;
	pDis = 0.;
	if (actBst.tgl != 0) {
		chkBst = FindElementArray(actBst.fil.l, actBst.fil.c, 
				act.iF[locActInd], 0, 3);
		if (chkBst > -1) {
			if (!(actBst.facKWA == 0. && cnt > 0) 
					&& P2A(actBst.fil.l,chkBst,1,3) == side) {
				pDis = (cnt > 0) ? actBst.facKWA : 1.;
			}
			else {chkBst = -1; }
		}
	}
	if (chkBst == -1 && actDis.tgl != 0) {
		pDis = (cnt > 0 ) ? actDis.pWA[side] : actDis.p[side];
	}

	pDis = AdjustDynamicsRate(pDis);
	CONT(!(genrand_real3() < pDis));
	if (act.fix[locActInd] > -1) { 
		act.fix[locActInd] = -1;
		if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) { act.nFA[locActInd] = 0; }
  		if (rheoWay > 0) { 
			DeleteElement1dArray(meaStrePar.l, &meaStrePar.c, actInd);
			DeleteElement1dArray(appStraPar.l, &appStraPar.c, actInd);
		}
	}
	ind = (side == 0) ? actInd2 : actInd;
	for(k = 2; k < nChAc; k++) {
		abpInd = P2A(act.ch,iAct[ind],k,nChAc);
		// If there is ABP bound on the actin
		CONT(abpInd < 0);
		UpdateActinDisassemblySubroutine2(abpInd, ind);
		V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, ind, k);
		(sendAbpDyn.c)++;
	} 

	// Change information of the disassembled particle
	UpdateActinAbpMonomerListSubroutine2(locActInd, actInd, 0);
	InsertElement1dArrayWoChk(actM.l, &actM.c, actInd);
	(actDis.cntMe)++;
	if (recTraj.tgl != 0) {
		ReplaceElementInTrajList(actInd, 0);
	}
	V3SET(&P2A(sendActDyn.l,sendActDyn.c,0,3), actInd, actInd2,
			-1 * (side + 1));
	(sendActDyn.c)++;

	UpdateActinDisassemblySubroutine(actInd2, side);
	// Update the neighboring list
	// Find and delete an actin cylinder corresponding the disassembled one
	DeleteElementInNeighborList(actInd, 0);

  	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), actInd, currTimeStep);
  	(noActDyn.c)++;
	n--;
  }
}

// Choose filaments which will undergo the bursting depolymerization.
void UpdateActinBurstDisassembly(void) {
}

void UpdateActinDegradation(int *ind, double len, int dim) {
}

// Eliminate expired elements in noActDyn.l
void UpdateNoActinDynamicsList(void) {
  int n;

  CheckArraySize(&noActDyn, &noActDyn.siz, 2, 1);
  for (n = noActDyn.c - 1; n >= 0; n--) {
    CONT(!(currTimeStep - P2A(noActDyn.l,n,1,2) > durNoActDyn));
	DeleteElementArrayByIndex(noActDyn.l,&noActDyn.c,n,2);
  }
}

int CheckActinAvailability(int actInd, int stage) {
  int CS;

  CS = FindElementArray(noActDyn.l, noActDyn.c, actInd, 0, 2);
  if (CS == -1 && stage >= 1) { 
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, actInd, 1, 3);
  }
  return CS;
}

/*------------------------ Dynamic behaviors of actins -----------------------*/

/*------------------------- Dynamic behaviors of ABPs ------------------------*/

//double UpdateMotorWalkingSubroutine(int abpInd, int actInd, int side) {
double UpdateMotorWalkingSubroutine(int abpInd, int side) {
  int fMag;
 
  fMag = motWalk.maxF[0] - 1
		+ (int)(P2A(recInstSprFabp,iAbp[abpInd],side * 2 + 1,4));
  // Assign walking probability depending on the longitudinal force
  fMag = TrimIntVal(fMag, 0, motWalk.maxF[0] + motWalk.maxF[1] - 2);
  return motWalk.p[fMag];
}

// Update the walking event of motor
void UpdateMotorWalking(void) {
  int n, k, l, loc, loc2, CS, *pArr, *pArr2, nextLoc, startLoc;
  int actInd, abpInd, locActInd, nextActInd, locNextActInd;
  double pMotW;
  ListInt confAct;

  confAct.l = allIntL;
  confAct.c = 0;

  FOR_ABPME(n) {
	abpInd = abp.id[n];
	pArr = &P2A(abp.ch,n,0,nChAb);
	// Only motors are involved in this function.
	CONT(pArr[2] != 2);
	// free motors are not considered in this function.
	CONT(pArr[0] < 0 && pArr[1] < 0);
	// If motor unbinds or binds at the current time step, it skips 
	// walking procedure
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
	CONT(CS > -1);

	CONT(abp.mId[n] < (int)(motSA.nMotPerSide * motSA.porInac));

    for (k = 0; k < 2; k++) {
		actInd = pArr[k];
		CONT(actInd < 0);
		CS = CheckActinAvailability(actInd, 1);
		CONT(CS > -1);
        locActInd = iAct[actInd];
		pArr2 = &P2A(act.ch,locActInd,0,nChAc);

		loc = FindAbpActinChain(locActInd, abpInd, 0);	
		loc2 = (loc - 2) / nChAcY;
		// If ABP is bound to the last set of binding sites, the walking 
		// motion makes ABP belong to the next actin segment in the barbed
		// direction
		if (loc2 == nChAcX - 1) {
			nextActInd = pArr2[0];
			locNextActInd = iAct[nextActInd];
			// If next binding site is a barbed end, motors cannot walk.
			CONT(P2A(act.ch,locNextActInd,0,nChAc) < 0);
			startLoc = 2;
			CS = CheckActinAvailability(nextActInd, 1);
			CONT(CS > -1);
		}
		else {
			nextActInd = actInd;
			locNextActInd = locActInd;
			startLoc = 2 + (loc2 + 1) * nChAcY;
		}
		for(l = 0; l < nChAcY; l++) {
			CONT(P2A(act.ch,locNextActInd,startLoc + l,nChAc) > -1);
			CS = Find2ElementArray(confAct.l, confAct.c, locNextActInd, 
					startLoc + l, 0, 2);
			CONT(CS > -1);
			nextLoc = startLoc + l;
			break;
		}
		// If there is not available binding site, motors cannot walk.
		CONT(l == nChAcY);
		// Find longitudinal force
		pMotW = UpdateMotorWalkingSubroutine(abpInd, k);
		pMotW = AdjustDynamicsRate(pMotW);
		// Generate random number to decide whether it walks or not
        CONT(!(genrand_real3() < pMotW));
		// Adjust chain of actin
		pArr2[loc] = -1;
		P2A(act.ch,locNextActInd,nextLoc,nChAc) = abpInd;
		if (actInd != nextActInd) {
			// Adjust chain of motor
			pArr[k] = nextActInd;
			if (mpiMethod == 0) {
				// Adjust longCh.l
				DeleteLongChain(abpInd + nAct, actInd);
				InsertLongChain(abpInd + nAct, nextActInd, minDimDomC * 0.9);
			}
		}
		(motWalk.cntMe)++;
		if (recAbpDyn.tgl != 0) { abpDyn.cntW[n]++; }
		// Put them in sendAbpDyn.l to let other CPUs know about this event
		V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, actInd, loc);
		V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c + 1,0,3), abpInd, 
				nextActInd, nextLoc);
		(sendAbpDyn.c) += 2;
		// Put Them in noAbpDyn.l
		V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), 
				abpInd, actInd, -1 * currTimeStep);
		V3SET(&P2A(noAbpDyn.l,noAbpDyn.c + 1,0,3), 
				abpInd, nextActInd, -1 * currTimeStep);
		(noAbpDyn.c) += 2;
		V2SET(&P2A(confAct.l,confAct.c,0,2), locActInd, loc);
		(confAct.c)++;
		break;
	}
  }
}
// Subroutine for UpdateAbpBinding()
void UpdateAbpBindingSubroutine(int *actInd, int abpInd) {
  int n, side, findSide, kind, CS;
  int *cntRebMe, *cntMoBindMe, *nAbpInaMe, *nAbpMme, *pArr, *pArr2;
  int locActInd[2], locAbpInd, tglRebAng, ind, locOppActInd;
  double dr1[NDIM], dr2[NDIM], drCrs[2][NDIM];
  double rPos[NDIM], rPos2[NDIM], rPos3[NDIM];
  double dp, len, pAbpBind, prevLen;

  locActInd[0] = iAct[actInd[0]];
  locActInd[1] = iAct[actInd[1]];
  
  locAbpInd = iAbp[abpInd];
  pArr = &P2A(act.ch,locActInd[0],0,nChAc);
  pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
  kind = K_ABP(locAbpInd);
  if (pArr2[0] > -1) {
	locOppActInd = iAct[pArr2[0]]; 
  }

  // Motor
  if (kind == 2) {
	pAbpBind = motReb.p;
	cntRebMe = &motReb.cntMe;
	cntMoBindMe = &motMoBind.cntMe;
	nAbpInaMe = &nMotInaMe;
	nAbpMme = &nMotMme;
	tglRebAng = motReb.gTglCrsAng;
  }
  // ACP
  else {
	pAbpBind = acpReb.p;
	cntRebMe = &acpReb.cntMe;
	cntMoBindMe = &acpMoBind.cntMe;
	nAbpInaMe = &nAcpInaMe;
	nAbpMme = &nAcpMme;
	tglRebAng = acpReb.gTglCrsAng;
  }
  pAbpBind = AdjustDynamicsRate(pAbpBind);
  CS = 1;
  // Prevent ABP from binding to the same filament twice
  if (pArr2[0] > -1) {
	if (act.iF[locActInd[0]] == act.iF[locOppActInd]) 
	{ CS = 0; }
  }
  // Check distance between center of ABP and actin
  if (CS == 1) {  
	CS = 0;
	CalcVec(dr1, &P2(act.r,locActInd[1],0), &P2(act.r,locActInd[0],0));
	V3SCALE(dr1, INV(nChAcX));
	V3COPY(rPos, &P2(act.r,locActInd[1],0));
	prevLen = POS_LARGE_VALUE;
	// Check available positions from the barbed direction because motors
	// walk toward the barbed end.
	for(n = nChAcX - 1; n >= 0; n--) {
		VV3SUB(rPos, dr1);
		ApplyBoundCondVector(rPos, -1, 0);	
		// One of ABP binding sites between two end points should be 
		// empty at least.
		findSide = FindElementArray(&pArr[n * nChAcY + 2], 
				nChAcY, -1, 0, 1); 
		CONT(findSide < 0);
		len = CalcDist(rPos, &P2(abp.r,locAbpInd,0), 0);
		if (ISMTF(kind)) {
			// Binding toward pointed ends is prevented.
			// Check the distance between binding point and ABP.
			CONT(!(len >= abpF.spr[kind].lo && len <= abpF.spr[kind].hi));
			// Find the nearest binding point from ABP.
			BREAK(prevLen < len);
			prevLen = len;
			side = findSide + n * nChAcY + 2;
			V3COPY(rPos2, rPos);
			CS = 1;
		}
		else {
			if (len >= abpF.spr[kind].lo && len <= abpF.spr[kind].hi) { 
				side = findSide + n * nChAcY + 2;
				CS = 1;
				break;
			}
		}
	}
  }
  if (!(ISMTF(kind))) { 
	// Check distance between the center of actin and the end point of ABP.
	// 0.1 is arbitrarily determined.
	if (CS == 1 && abpF.bend[kind].stf > 0. && pArr2[0] > -1) {
		CalcAbpArmEndPos(rPos3, pArr2[0], abpInd);
		len = CalcDist(rPos, rPos3, 0);
	    if (!(len < 0.5 * actF.dia + 0.1)) { CS = 0; }
	}
	// Check angle formed by axis of act.filament and potential ABP chain.
	// The angle should be closer to right angle
	if (CS == 1 && abpF.a90[kind].stf > 0.) {
		CalcVec(dr1, &P2(act.r,locActInd[1],0), &P2(act.r,locActInd[0],0));
		CalcVec(dr2, rPos, &P2(abp.r,locAbpInd,0));
	    dp = fabs(V3COS(dr1, dr2));
	    if (!(dp < abpF.a90[kind].lo)) { CS = 0; }
	}
  }
  else {
	// Prevent the binding of binding sites on the same side.
	// This is for reflecting the structure of myosin thick filaments where
	// two heads from one myosin face in the opposite directions.
	if (CS == 1 && pArr2[0] > -1 && motReb.gTglOppDir != 0) {
		if (locOppActInd > -1) {
			CalcVec(dr1, &P2(act.r,locActInd[1],0), &P2(act.r,locActInd[0],0));
			CalcVecActinAbp(dr2, pArr2[0], abpInd, 0);
			V3CROSS(drCrs[0], dr1, dr2);
			CalcVec(dr2, rPos2, &P2(abp.r,locAbpInd,0));
			V3CROSS(drCrs[1], dr1, dr2);
			if (V3DOT(drCrs[0], drCrs[1]) > 0) {
				CS = 0;
			}
		}
	}
	// Allow the binding only when thick filaments are aligned properly along
	// with the polarity of actin filaments.
	if (CS == 1 && motReb.gTglCrsAng != 0) {
		CS = 0;
		if (pArr2[3] > -1 || pArr2[4] > -1) {
			ind = (pArr2[3] > -1) ? iAbp[pArr2[3]] : iAbp[pArr2[4]];
			CalcVec(dr1, &P2(act.r,locActInd[1],0), &P2(act.r,locActInd[0],0));
			CalcVec(dr2, &P2(abp.r,locAbpInd,0), &P2(abp.r,ind,0));
			if (abp.mId[locAbpInd] - abp.mId[ind] < 0) {
				V3REVSIGN(dr2);
			}
			if (V3DOT(dr1, dr2) > 0) {
				CS = 1;
			}
		}
	}
  }
  // Check angle formed by the axes of two actin filaments
  if (CS == 1 && (!(ISMTF(pArr2[2]))) && tglRebAng != 0 && pArr2[0] > -1) {
	if (locOppActInd > -1) {
		CalcVec(dr1,&P2(act.r,locActInd[1],0), &P2(act.r,iAct[pArr[1]],0));
		V3SUB(dr2, &P2(act.r,iAct[P2A(act.ch,locOppActInd,0,nChAc)],0),
				&P2(act.r,iAct[P2A(act.ch,locOppActInd,1,nChAc)],0));
		ApplyBoundCondVecDiff(dr2);
		dp = V3COS(dr1, dr2);
		if (kind == 1 && dp < 0) { dp = REVSIGN(dp); }
		if (kind == 2 && dp > 0 && tglRebAng == 2) 
		{ dp = REVSIGN(dp); }
		if (dp < abpF.cr[kind].lo || dp > abpF.cr[kind].hi) { CS = 0; } 
	}
	else { CS = 0; }
  }
 
  if (CS == 1 && genrand_real3() < pAbpBind) {
	pArr[side] = abpInd;
	if (recAbpDyn.tgl != 0) { abpDyn.cntB[locAbpInd]++; }
	if (pArr2[0] < 0 && pArr2[1] < 0) { 
		(*cntMoBindMe)++;
		pArr2[0] = actInd[0];
		(*nAbpMme)--;
		(*nAbpInaMe)++;
	}
	else {
		(*cntRebMe)++;
		pArr2[1] = actInd[0];
		(*nAbpInaMe)--;
		if (tglNeiAbpSC != 0 && tglNeiAbpDC == 0) {
			DeleteElementInNeighborList(abpInd, 1);
		}
	}
	UpdateAbpUnbRebLists(abpInd, actInd[0], side, 1);
	if (mpiMethod == 0) {
		InsertLongChain(abpInd + nAct, actInd[0], minDimDomC * 0.9);
	}
  }
}
// Update the binding event of ACP or motor
void UpdateAbpBinding(void) {
  int n, CS, *nlPnt, abpInd, locAbpInd, *actCyl;

  nlPnt = neigh.l;
  for(n = 0; n < neigh.c; n++) {
    if (n > 0) { nlPnt += 2; }
	// Consider the possible reformation of ACP chains
	CONT(!((nlPnt[0] >= nAct && nlPnt[1] < nAct) 
			|| (nlPnt[0] < nAct && nlPnt[1] >= nAct)));
	actCyl = &P2A(act.cyl.l,nlPnt[(nlPnt[0] < nAct) ? 0 : 1],0,2);
	abpInd = nlPnt[(nlPnt[0] < nAct) ? 1 : 0] - nAct;
	locAbpInd = iAbp[abpInd];
	CONT(locAbpInd < 0 || locAbpInd >= nAbpMe);
	CONT(iAct[actCyl[0]] < 0 || iAct[actCyl[1]] < 0);
	// One of the actin binding site should be available.
	CONT(P2A(abp.ch,locAbpInd,0,nChAb) > -1 
			&& P2A(abp.ch,locAbpInd,1,nChAb) > -1);
	// Selective control of binding depending on the kind of ABPs.
	CONT((K_ABP(locAbpInd) == 2 && motReb.tgl == 0) 
			|| (K_ABP(locAbpInd) != 2 && acpReb.tgl == 0));
	CS = CheckActinAvailability(actCyl[0], 1);
	CONT(CS > -1);
	// Prevent the unbound ABP from reforming during a certain time 
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
	CONT(CS > -1);
	UpdateAbpBindingSubroutine(actCyl, abpInd);
  }
}

void UpdateAbpMonomerBinding(void) {
  int n, k, *pArr, side, ind, CS, kind, rankMol, cntTry; 
  int abpInd, locAbpInd;
  double dimDomC[NDIM], rPos[NDIM], rOri[NDIM];
  double pAcpMoBind, pMotMoBind, randNum, vol;
  ListInt *pM;

  V3SUB(dimDomC, &P2A(bnd.r,1,0,NDIM), &P2A(bnd.r,0,0,NDIM));
  vol = V3PROD(dimDomC);
  if (gTglImpAcpM != 0) {
	  pAcpMoBind = 1. - exp(acpMoBind.facP * (double)acpM.c / vol);
	  pAcpMoBind = AdjustDynamicsRate(pAcpMoBind) / (nChAc - 2);
  }
  else { pAcpMoBind = 0.; }
  if (gTglImpMotM != 0) {
	pMotMoBind = 1. - exp(motMoBind.facP * (double)motM.c / vol);
	pMotMoBind = AdjustDynamicsRate(pMotMoBind) / (nChAc - 2);
  }
  else { pMotMoBind = 0.; }

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	BREAK(!(pAcpMoBind + pMotMoBind > 0.));
	if (acpM.c == 0) { pAcpMoBind = 0.; }
	if (motM.c == 0 || (motSA.gTgl != 0 && motSA.nNucMe == 0)) { 
		pMotMoBind = 0.;
	}
	// Check whether it is ready for dynamics.
	CS = CheckActinAvailability(act.id[n], 1);
	CONT(CS > -1);
	pArr = &P2A(act.ch,n,0,nChAc);
	// If it is a barbed end, the binding cannot occur.
	CONT(pArr[0] < 0);
	// Choose a random position from which available site will be searched.
	side = GenRandIntIndex(nChAc - 2) + 2;
	for(k = 0; k < nChAc - 2; k++) {
		BREAK(P2A(act.ch,n,side,nChAc) < 0);
		side++;
		if (side == nChAc) { side = 2; }
	}
	CONT(k == nChAc - 2);

	// Check a probability
	randNum = genrand_real3();
	CONT(!(randNum < pAcpMoBind + pMotMoBind));
	// Find a monomeric ABP in the list
	pM = (randNum < pMotMoBind) ? &motM : &acpM;
	for(ind = 0; ind < pM->c; ind++) {
		abpInd = pM->l[ind];
		// Check whether the motor is ready for dynamics..
		CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
		BREAK(CS == -1);
	}
	//BREAK(ind == pM->c);
	CONT(ind == pM->c);
	locAbpInd = iAbp[abpInd];
	kind = K_ABP(locAbpInd);
	// Calculate the new position of the ABP
	cntTry = 0;
	while(1) {
		cntTry++;
		// In a certain location, it might not be possible to accommodate 
		// an ABP. An infinite number of trials should be prevented.
        BREAK(cntTry == 1000);
		CalcInactAbpPosition(rPos, rOri, &P2(act.r,n,0),
				&P2(act.r,iAct[pArr[0]],0), kind, side);
		if (bnd.gTglRnd != 0) {
			CS = CheckActinAbpOverlapBoundary(rPos);
			CONT(CS != 1);
		}
		CS = CheckParticleInDomain(rPos);
		CONT(CS != 1);
		break;
    }
    CONT(cntTry == 1000);
	DeleteElementArrayByIndex(pM->l, &pM->c, ind, 1);
	if (kind == 2) { 
		pMotMoBind = 1. - exp(motMoBind.facP * (double)motM.c / vol);
		pMotMoBind = AdjustDynamicsRate(pMotMoBind) / (nChAc - 2);
	}
	else {
		pAcpMoBind = 1. - exp(acpMoBind.facP * (double)acpM.c / vol);
		pAcpMoBind = AdjustDynamicsRate(pAcpMoBind) / (nChAc - 2);
	}
	V3COPY(&P2(abp.r,locAbpInd,0), rPos);
	V3SET_ALL(&P2(abp.fBr,locAbpInd,0), 0.);
	// Update the chain information of an ABP
	pArr[side] = abpInd;
	P2A(abp.ch,locAbpInd,0,nChAb) = act.id[n];
	if (abpAge.gTgl != 0) { abp.age[locAbpInd] = currTimeStep; }
	// Adjust the counter to reflect the binding
	if (kind == 2) {
		if (motSA.gTgl != 0) { 
			abp.mId[locAbpInd] = 0;
			motSA.cntNucMe++; 
			motSA.nNucMe--;
		}
		nMotMme--;
		nMotInaMe++;
		(motMoBind.cntMe)++;
	}
	else {
		nAcpMme--;
		nAcpInaMe++;
		(acpMoBind.cntMe)++;
	}
	if (recAbpDyn.tgl != 0) { abpDyn.cntB[locAbpInd]++; }
	UpdateAbpUnbRebLists(abpInd, act.id[n], side, 0);
	if (mpiMethod == 0) {
		InsertLongChain(abpInd + nAct, act.id[n], minDimDomC * 0.9);
	}
	if (tglNeiAbpSC != 0) {
		InsertElementInNeighborList(abpInd, -1, 2);
	}
  }
}

void UpdateMotorTurnover(void) {
  int n, k, CS, curr, next, locAbpInd, actInd, locActInd, actSide, cntAct;
  int *pArr, *pArr2, fMag, cntF;
  double pMotTurn, len, f, lenEq;

  FOR_ABPME(n) { 
	pArr = &P2A(abp.ch,n,0,nChAb);
	// It should be a motor
	CONT(pArr[2] != 2);
	// It should be the end of motor filaments
	CONT(!(pArr[3] > -1 && pArr[4] < 0));
	
	cntF = 0;
	f = 0.;
	curr = abp.id[n];
	CS = 1;
	cntAct = 0;
	while(curr > -1) {
		locAbpInd = iAbp[curr];
		if (locAbpInd < 0 || locAbpInd >= nAbpMe) { 
			CS = 0;
			break;
		}
		pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
		next = pArr2[3];
		if (next > -1 && motSA.to.mode == 1) {
	        len = CalcDist(&P2(abp.r,locAbpInd,0), &P2(abp.r,iAbp[next],0), 0);
	        lenEq = (abp.mId[locAbpInd] == 0 && abp.mId[iAbp[next]] == 0)
	                ? motSA.cenDist : motSA.spr.eq;
	        f += motSA.spr.stf * (len - lenEq);
			cntF++;
		}
		for(k = 0; k < 2; k++) {
			actInd = pArr2[k];
			CONT(actInd < 0);
			if (motSA.to.mode == 0) {
				CS = 0;
				break;
			}
			else {
				cntAct++;
				locActInd = iAct[actInd];
				if (locActInd < 0 || locActInd >= nActMe) {
					CS = 0;
					break;
				}
			}
		}
		BREAK(CS == 0);
		curr = next;
	}
	CONT(CS != 1);
	
	if (motSA.to.mode == 1) { 
		if (cntAct == 0) {
			pMotTurn = motSA.to.p[0];
		}
		else {
			if (cntF > 0) { f /= (double)cntF; }
			if (f < 0.) { f = 0.; }
			fMag = TrimIntVal((int)f, 0, motSA.to.maxF - 1);
			pMotTurn = motSA.to.p[fMag];
		}
	}
	else {
		pMotTurn = motSA.to.pf;
	}
	CONT(!(genrand_real3() < pMotTurn));
	
	curr = abp.id[n];
	while(curr > -1) {
		locAbpInd = iAbp[curr];
		pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
		next = pArr2[3];
		if (tglNeiAbpSC != 0 || tglNeiAbpDC != 0) {
			DeleteElementInNeighborList(curr, 1);
		}
		cntAct = 0;
		for(k = 0; k < 2; k++) {
			actInd = pArr2[k];
			CONT(actInd < 0);
			cntAct++;
			locActInd = iAct[actInd];
			// UNBINDING			
  			actSide = FindAbpActinChain(locActInd, curr, 0);
			P2A(act.ch,locActInd,actSide,nChAc) = -1;
			if (mpiMethod == 0) { 
				DeleteLongChain(curr + nAct, actInd);
			}
			V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), curr, actInd, 
					currTimeStep);
			(noAbpDyn.c)++;
	    	(motUnb.cntMe)++;
			pArr2[k] = -1; 
		}
		if (mpiMethod == 0 && next > -1 && motSA.cenDist * 1.5 > neiEdge) {
			if (abp.mId[locAbpInd] == 0 && abp.mId[iAbp[next]] == 0) {
				DeleteLongChain(curr + nAct, next + nAct);
			}
		}
		V2SET_ALL(&pArr2[3], -1);
		V3SET_ALL(&P2(abp.r,locAbpInd,0), 0.);
		V3SET_ALL(&P2(abp.f,locAbpInd,0), 0.);
		V3SET_ALL(&P2(abp.fBr,locAbpInd,0), 0.);
		abp.mId[locAbpInd] = -1;
		if (abpAge.gTgl != 0) { abp.age[locAbpInd] = 0; }
		InsertElement1dArrayWoChk(motM.l, &motM.c, curr);
		if (recTraj.tgl2 != 0) {
			ReplaceElementInTrajList(curr, 1);
		}
		if (cntAct == 0) {
			V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), curr, -1, currTimeStep);
			(noAbpDyn.c)++;
		}
		else {
	    	(nMotMme)++;
		}
		if (cntAct == 1) {
			nMotInaMe--;
		}
		if (recAbpDyn.tgl != 0) { abpDyn.cntU[locAbpInd]++; }
		curr = next;
	}
	motSA.cntTurnMe++;
	(motSA.nNucMe)++;
  }
}

// Calculate self-assembly of motors 
void UpdateMotorAssembly(void) {
  int n, side, ind, cntTry, CS, *pArr, abpInd, locAbpInd;
  double pMotAss, dr[NDIM], rPos[NDIM];

  pMotAss = AdjustDynamicsRate(motSA.pAss);
  FOR_ABPME(n) {
	CONT(ISABPM(n));
	BREAK(motM.c == 0);
	pArr = &P2A(abp.ch,n,0,nChAb);
	// It should be a motor
	CONT(pArr[2] != 2);
	// It should be one of the ends of motor filaments
	CONT(pArr[3] > -1 && pArr[4] > -1);
	// Check when the previous dynamic event occurred for motor
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abp.id[n], 0, 3);
	CONT(CS > -1);
	// Check a probability
	CONT(!(genrand_real3() < pMotAss));
	// Determine the side for assembly
	side = (pArr[3] < 0) ? 0 : 1;
	// Find a monomeric motor
	for(ind = 0; ind < motM.c; ind++) {
		abpInd = motM.l[ind];
		CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
		BREAK(CS == -1);
	}
	BREAK(ind == motM.c);
	locAbpInd = iAbp[abpInd];

	// If it is a seed, a bare zone with large spacing is formed first.
	if (pArr[4 - side] < 0) {
		CS = 0;
		cntTry = 0;
		// If it is out of boundary, nucleation cannot happen
		while (CS != 1) {
			cntTry++;
	        BREAK(cntTry == 1000);
			GenRandDirecVec(dr);
			VS3ADD(rPos, &P2(abp.r,n,0), dr, motSA.cenDist);
			ApplyBoundCondVector(rPos, -1, 0);
			CS = CheckParticleInDomain(rPos);
			CONT(CS != 1);
			if (bnd.gTglRnd != 0) {
				CS = CheckActinAbpOverlapBoundary(rPos);
				CONT(CS != 1);
			}
		}
		CONT(cntTry == 1000);
	}
	else {
		CalcUnitVec(dr,&P2(abp.r,n,0),
				&P2(abp.r,iAbp[pArr[4 - side]],0));
		VS3ADD(rPos, &P2(abp.r,n,0), dr, motSA.spr.eq);
		ApplyBoundCondVector(rPos, -1, 0);
		CS = CheckParticleInDomain(rPos);
		CONT(CS != 1);
		if (bnd.gTglRnd != 0) {
			CS = CheckActinAbpOverlapBoundary(rPos);
			CONT(CS != 1);
		}
		// If the size of motor filaments should be the same to each other
		if (motSA.gTglConSiz != 0) {
			CONT(abp.mId[n] >= motSA.nMotPerSide - 1);
		}
	}
	DeleteElementArrayByIndex(motM.l, &motM.c, ind, 1);
	V3COPY(&P2(abp.r,locAbpInd,0),rPos);
	V3SET_ALL(&P2(abp.fBr,locAbpInd,0), 0.);
	P2A(abp.ch,n,side + 3,nChAb) = abpInd;
	P2A(abp.ch,locAbpInd,4 - side,nChAb) = abp.id[n];

	motSA.cntAssMe++;
	if (abpAge.gTgl != 0) { abp.age[locAbpInd] = currTimeStep; }
	if (mpiMethod == 0 && pArr[4 - side] < 0 && motSA.cenDist * 1.5 > neiEdge) {
		InsertLongChain(abp.id[n] + nAct, abpInd + nAct, minDimDomC * 0.9);
	}
	abp.mId[locAbpInd] = (pArr[4 - side] < 0) ? 0 
			: abp.mId[n] + 1;
	V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abp.id[n], -1, currTimeStep);
	V3SET(&P2A(noAbpDyn.l,noAbpDyn.c + 1,0,3), abpInd, -1, currTimeStep);
	(noAbpDyn.c) += 2;
	if (tglNeiAbpSC != 0) {
		InsertElementInNeighborList(abpInd, -1, 2);
	}
  }
}

double CalcUnbindingRate(int abpInd, int actInd, int abpSide, int actSide) {
  int fMag, locActInd, locAbpInd;
  double pUnb;

  locActInd = iAct[actInd];
  locAbpInd = iAbp[abpInd];
  // If it is a motor, use the longitudinal force.
  if (K_ABP(locAbpInd) == 2) {
	fMag = motUnb.maxF[0] - 1 + (int)(P2A(recInstSprFabp,locAbpInd,
			abpSide * 2 + 1,4));
	fMag = TrimIntVal(fMag, 0, motUnb.maxF[0] + motUnb.maxF[1] - 2);
	pUnb = motUnb.p[fMag];
	if (gTglMotWalkSld != 0 && motWalk.tgl != 0) {
		if (P2A(act.ch,iAct[P2A(act.ch,locActInd,0,nChAc)],0,nChAc) < 0
				&& (int)((actSide - 2) / nChAcY) == nChAcX - 1) {
			fMag = motWalk.maxF[0] - 1
					+ (int)P2A(recInstSprFabp,locAbpInd,abpSide * 2 + 1,4);
			fMag = TrimIntVal(fMag, 0, motWalk.maxF[0] + motWalk.maxF[1] - 2);
			pUnb += motWalk.p[fMag];
		}
	}
  }
  // If it is an ACP, just use a force.
  else {
	fMag = (int)P2A(recInstSprFabp,locAbpInd,abpSide * 2,4);
	fMag = TrimIntVal(fMag, 0, acpUnb.maxF - 1);
	pUnb = acpUnb.p[fMag];
  }
  pUnb = AdjustDynamicsRate(pUnb);

  return pUnb;
}

void UpdateActiveAbpUnbindingSubroutine(int actInd, int abpInd) {
  int k, locAbpInd, locActInd, locOppActInd, abpSide, actSide;
  int *nAbpInaMe, *cntUnbMe, *pArr;

  locAbpInd = iAbp[abpInd];
  locActInd = iAct[actInd];
  pArr = &P2A(abp.ch,locAbpInd,0,nChAb);
  abpSide = (pArr[0] == actInd) ? 0 : 1;

  // Adjust chain of actin
  actSide = FindAbpActinChain(locActInd, abpInd, 0);
  P2A(act.ch,locActInd,actSide,nChAc) = -1;
  // Adjust chain of ABP
  if (abpSide == 0) { pArr[0] = pArr[1]; }
  pArr[1] = -1;

  if (locAbpInd < nAbpMe) { 
	nAbpInaMe = (K_ABP(locAbpInd) == 2) ? &nMotInaMe : &nAcpInaMe;
	cntUnbMe = (K_ABP(locAbpInd) == 2) ? &motUnb.cntMe : &acpUnb.cntMe;
    (*nAbpInaMe)++;
    (*cntUnbMe)++;
	if (recAbpDyn.tgl != 0) { abpDyn.cntU[locAbpInd]++; }

  }
  // The pair should be deleted in longCh.l
  if (locActInd < nActMe || locAbpInd < nAbpMe) {
	// Prevent the unbound ABP from reforming during a certain time.
	V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, currTimeStep);
	(noAbpDyn.c)++;
	if (mpiMethod == 0) { 
		DeleteLongChain(abpInd + nAct, actInd);
	}
  }
}

// Update the unbinding event of ACP or motor
void UpdateActiveAbpUnbinding(void) {
  int n, k, CS, actInd, locActInd, abpInd, side, *pArr;
  double pUnb;

  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
    abpInd = abp.id[n];
	CONT((K_ABP(n) == 2 && motUnb.tgl == 0) 
			|| (K_ABP(n) != 2 && acpUnb.tgl == 0));
	// Check noAbpDyn.l
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
	CONT(CS > -1);
    for (k = 0; k < 2; k++) {
		actInd = pArr[k];
		locActInd = iAct[actInd];
		CS = CheckActinAvailability(actInd, 1);
		CONT(CS > -1);
		side = FindAbpActinChain(locActInd, abpInd, 0);
		pUnb = CalcUnbindingRate(abpInd, actInd, k, side);
        CONT(!(genrand_real3() < pUnb));
		UpdateActiveAbpUnbindingSubroutine(actInd, abpInd);
		if (tglNeiAbpSC != 0 && tglNeiAbpDC == 0) {
			InsertElementInNeighborList(abpInd, -1, 2);
		}

		V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, actInd, side);
		(sendAbpDyn.c)++;
		n--;
		break;
	}
  }
}
void UpdateInactAbpUnbindingSubroutine(int actInd, int abpInd) {
  int locActInd, locAbpInd, side, kind, *pArr;
  ListInt *pM;
  
  locActInd = iAct[actInd];
  locAbpInd = iAbp[abpInd];
  pArr = &P2A(abp.ch,locAbpInd,0,nChAb);
  kind = K_ABP(locAbpInd);
  side = FindAbpActinChain(locActInd, abpInd, 0);
  P2A(act.ch,locActInd,side,nChAc) = -1;
  pArr[0] = -1;
  if ((kind != 2 && gTglImpAcpM != 0) 
		|| (kind == 2 && gTglImpMotM != 0 && motSA.gTgl == 0)) {
	V3SET_ALL(&P2(abp.r,locAbpInd,0), 0.);
  }
  if (locAbpInd < nAbpMe) {
	if ((kind != 2 && gTglImpAcpM != 0) 
			|| (kind == 2 && gTglImpMotM != 0 && motSA.gTgl == 0)) {
		pM = (kind == 2) ? &motM : &acpM;
		InsertElement1dArrayWoChk(pM->l, &pM->c, abpInd);
		if (recTraj.tgl2 != 0) {
			ReplaceElementInTrajList(abpInd, 1);
		}
	}
	if (!(ISMTF(kind) && (pArr[3] > -1 || pArr[4] > -1)) && abpAge.gTgl != 0) {
		abp.age[locAbpInd] = 0;
	}
	// Adjust counters
	// If motor
	if (kind == 2) {
		nMotInaMe--;
		nMotMme++;
		(motInaUnb.cntMe)++;
		if (motSA.gTgl != 0 && !(pArr[3] > -1 || pArr[4] > -1)) {
			motSA.nNucMe++;
		}
	}	
	// If ACP
	else {
		nAcpInaMe--;
		nAcpMme++;
		(acpInaUnb.cntMe)++;
	}
	if (recAbpDyn.tgl != 0) { abpDyn.cntU[locAbpInd]++; }

  }
  // The pair should be deleted in longCh.l
  if (locActInd < nActMe || locAbpInd < nAbpMe) {
	V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, currTimeStep);
	(noAbpDyn.c)++;
	if (mpiMethod == 0) { 
		DeleteLongChain(abpInd + nAct, actInd);
	}
  }
}

// Update the unbinding event of ACP or motor
void UpdateInactAbpUnbinding(void) {
  int n, CS, side, actInd, locActInd, abpInd, kind, *pArr;
  double pUnb;

  FOR_ABPME(n) {
    pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] < 0));
    abpInd = abp.id[n];
	kind = K_ABP(n);
	CONT((kind == 2 && motInaUnb.tgl == 0) 
			|| (kind != 2 && acpInaUnb.tgl == 0));
	actInd = pArr[0];
	CS = CheckActinAvailability(actInd, 1);
	CONT(CS > -1);
	// If motors can multimerize, they cannot go back to a monomer state..
	if (ISMTF(kind)) {
		CONT(pArr[3] < 0 && pArr[4] < 0);
	}
	// Check noAbpDyn.l
	CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, abpInd, 0, 3);
	CONT(CS > -1);
	locActInd = iAct[actInd];
	side = FindAbpActinChain(locActInd, abpInd, 0);
	if (kind == 2) {
		// If it is the arm of a multimerized motor, tension can exist.
		// Also, the sliding-off at barbed ends should be considered.
		pUnb = CalcUnbindingRate(abpInd, actInd, 0, side);
	}
	else {
		// If it is not the arm of a multimerized motor, there is no tension on
		// the arm.
		pUnb = AdjustDynamicsRate(acpUnb.p[0]);
	}

	CONT(!(genrand_real3() < pUnb));
	// Delete the element in the neighboring list
	if (tglNeiAbpSC != 0) {
		if ((kind != 2 && gTglImpAcpM != 0) 
				|| (kind == 2 && gTglImpMotM != 0 && motSA.gTgl == 0)
				|| (kind == 2 && motSA.gTgl != 0 && pArr[3] < 0 
				&& pArr[4] < 0)) {
			DeleteElementInNeighborList(abpInd, 1);
		}
	}
	UpdateInactAbpUnbindingSubroutine(actInd, abpInd);

	V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, actInd, side);
	(sendAbpDyn.c)++;
	n--;
  }
}

// Eliminate expired elements in noAbpDyn.l
void UpdateNoAbpUnbindList(void) {
  int n, ind, dur, timeStep;

  CheckArraySize(&noAbpDyn, &noAbpDyn.siz, 3, 1);
  for (n = noAbpDyn.c - 1; n >= 0; n--) {
	ind = iAbp[P2A(noAbpDyn.l,n,0,3)];
	timeStep = P2A(noAbpDyn.l,n,2,3);
	if (ind > -1 && ind < nAbpMe) {
		if (K_ABP(ind) != 2) { dur = durNoAcpUnbReb; }
		else {
			dur = (timeStep > 0) ? durNoMotUnbReb : durNoMotWalk;
			if (timeStep < 0) { timeStep = -1 * timeStep; }
		}
	}
	else { dur = 0; }
    CONT(!(currTimeStep - timeStep > dur));
	DeleteElementArrayByIndex(noAbpDyn.l, &noAbpDyn.c, n, 3);
  }
}

/*------------------------- Dynamic behaviors of ABPs ------------------------*/

/*-------------------------- Updating chains and lists -----------------------*/

void UpdateChainList(void) {
  int n, *pArr;

  nAcpInaMe = 0;
  nMotInaMe = 0;
  nAcpMme = 0;
  nMotMme = 0;
  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
    if (pArr[0] > -1 && pArr[1] < 0) {
		if (K_ABP(n) == 2) { nMotInaMe++; }
		else { nAcpInaMe++; }
    }
	else if (pArr[0] < 0 && pArr[1] < 0) {
		if (K_ABP(n) == 2) { nMotMme++; }
		else { nAcpMme++; }
	}
  }
}

/*-------------------------- Updating chains and lists -----------------------*/


void CalcBoundSpringForces(void) {
}

void CalcBoundRepulsiveForces(void) {
}

void UpdateBoundaryActinUnbindMature(void) {
}

void UpdateBoundaryActinBinding(void) {
}

void UpdateBoundaryLocation(void) {
}

/*------------------------ Handle boundary conditions ------------------------*/

// Check whether the chain crosses a boundary or not.
void CheckCrossBound(int *sft, double *dr, double *dr2) {
  int k;
 
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	CONT(!(fabs(dr2[k] - dr[k]) > 1.5 * dimDomH[k]));
	sft[k] += (dr2[k] > dr[k]) ? -1 : 1;
  }
}

// Apply periodic boundary condition (PBC) to a chain vector.
// If the absolute value of a component of the vector is greater than half of 
// domain width with PBC, it is added or subtracted by domain width.
void ApplyBoundCondVecDiff(double *dr) {
  int k;
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	if (dr[k] >= dimDomH[k]) { dr[k] -= dimDom[k]; }
	else if (dr[k] <= REVSIGN(dimDomH[k])) { dr[k] += dimDom[k]; }
  }
}

// If mode is -1, neglect repulsive condition.
// Otherwise, if mode is 0, actin. If mode is 1-3, ABP.
// If 5, nucleus. If 4, membrane
void ApplyBoundCondVector(double *r, int mode, int mode2) {
  int k, indBndF;
  double disp, dispAll, repF, drag, sft;

  FOR_NDIM(k) {
	disp = (k == dirStr && rheoWay > 0 && bulkRheoType == 0) 
			? r[dirNoPBC] * stra.acc : 0.;
	// If PBC
    if (pbc[k] == 1) {
		if (r[k] < disp + rGrid[k][0])
		{ r[k] += dimDom[k]; }
		else if (r[k] >= disp + rGrid[k][nGrid[k] - 1])
		{ r[k] -= dimDom[k]; }
	}
	// If no PBC and mode >= 0, repulsive condition is applied.
	else if (pbc[k] != 1 && mode >= 0) {
		sft = 0.;
		if (mode == 0) { drag = 1.; }
		else if (mode >= 1 && mode <= 3) { drag = abpF.drag[mode - 1].inv; }
		else if (mode == 4) { drag = INV(memb.dragR); }
		else if (mode == 5) { drag = INV(memb.nucDragR); }
		else { 
			drag = INV(bead.drag[mode - 6]); 
			sft = bead.rad[mode - 6];
		}
		if ((r[k] < disp + rGrid[k][0] + sft) 
				|| (r[k] >= disp + rGrid[k][nGrid[k] - 1] - sft)) {
			if (r[k] < disp + rGrid[k][0] + sft) {
				dispAll = disp + rGrid[k][0] + sft;
				indBndF = 0;
			}
			else {
				dispAll = disp + rGrid[k][nGrid[k] - 1] - sft;
				indBndF = 1;
			}
			repF = bnd.stfRep * (r[k] - dispAll);
			// If shear-bulk-rheology (against a tilted boundary)
			if (k == dirStr && rheoWay > 0 && bulkRheoType == 0) {
				r[dirStr] -= repF / (SQR(stra.acc) + 1.) * drag * dt;
				r[dirNoPBC] -= (repF  * stra.acc / (SQR(stra.acc) + 1.)) 
						* drag * dt;
				if (mode2 != 0) {
					P2A(bnd.f,2 * dirStr + indBndF,dirStr,NDIM) += repF 
							/ (SQR(stra.acc) + 1.);
					P2A(bnd.f,2 * dirNoPBC + indBndF,dirNoPBC,NDIM) 
							+= repF * stra.acc / (SQR(stra.acc) + 1.);
				}
			}
			// If a flat boundary
			else {
				r[k] -= repF * drag * dt;
				if (bndMv.gTgl != 0 && mode2 != 0) {
					P2A(bnd.f,2 * k + indBndF,k,NDIM) += repF;
				}
			}
		}
	}
  }
}

// Apply the periodic boundary condition on the positions of whole things
void ApplyBoundCondAll(void) {
  int n;

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	ApplyBoundCondVector(&P2(act.r,n,0), 0, (act.fix[n] > -1 ? 0 : 1));
  }
 
  FOR_ABPME(n) {
	CONT(ISABPIM(n));
	ApplyBoundCondVector(&P2(abp.r,n,0), K_ABP(n) + 1, 1);
  }
}

// Regardless of shear deformation, offset positions of all particles into 
// retangular-solid domain
// If mode = 0, adjust the position only in shearing direction
// If mode > 0, apply PBC in other directions.
void ConvertRectDomainVector(double *r, int mode) {
  int k;

  if (mode == 0) {
	if (rheoWay > 0) {
		if (r[dirStr] < rGrid[dirStr][0]) { 
			while(r[dirStr] < rGrid[dirStr][0]) { 
				r[dirStr] += dimDom[dirStr];
			}
		}
		else if (r[dirStr] >= rGrid[dirStr][nGrid[dirStr] - 1]) {
			while(r[dirStr] >= rGrid[dirStr][nGrid[dirStr] - 1]) {
				r[dirStr] -= dimDom[dirStr];
			}
		}
	}
  }
  else {
	FOR_NDIM(k) {
		if (pbc[k] != 0 || (pbc[k] == 0 && confPbc[k] != 0)) {
			if (r[k] >= rGrid[k][nGrid[k] - 1]) { r[k] -= dimDom[k]; }
			else if (r[k] < rGrid[k][0]) { r[k] += dimDom[k]; }
		}
	}
  }
}

/*------------------------ Handle boundary conditions ------------------------*/

/*---------------- Interactions between boundaries and others ----------------*/

int CheckParticleInDomain(double *r) {
  int CS, k;
  double disp;
  CS = 1;
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 0));
	disp = (k == dirStr && rheoWay > 0 && bulkRheoType == 0) 
			? r[dirNoPBC] * stra.acc : 0.; 
	CONT(!(r[k] < disp + rGrid[k][0] || r[k] >= disp + rGrid[k][nGrid[k] - 1]));
	CS = 0;
	break;
  }
  return CS;
}

double HowManyInOutBound(void) {
  int k, cnt, dirX, dirY;
  int ind[NDIM], begin[NDIM], end[NDIM];
  double dist, resol, r[NDIM], cenDom[NDIM];

  resol = 0.5;
  FOR_NDIM(k) {
    cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
  }

  FOR_NDIM(k) {
    begin[k] = (int)ceil(P2A(bnd.r,0,k,NDIM) / resol);
    end[k] = (int)floor(P2A(bnd.r,1,k,NDIM) / resol);
  }

  if (dir2D == -1) {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			for(ind[2] = begin[2]; ind[2] < end[2]; ind[2]++) {
				VS3COPY(r, ind, resol);
				dist = CalcDist(cenDom, r, 0);
				if (dist <= bnd.radRnd) {
					cnt++;
				}
			}
		}
	}
	return (int)(cnt * CUBE(resol));
  }
  else {
	dirX = (dir2D + 1) % NDIM;
	dirY = (dir2D + 2) % NDIM;
	cnt = 0;
	for(ind[dirX] = begin[dirX]; ind[dirX] < end[dirX]; ind[dirX]++) {
		for(ind[dirY] = begin[dirY]; ind[dirY] < end[dirY]; ind[dirY]++) {
			r[dirX] = ind[dirX] * resol;
			r[dirY] = ind[dirY] * resol;
			r[dir2D] = cenDom[dir2D];
			dist = CalcDist(cenDom, r, 0);
			if (dist <= bnd.radRnd) {
				cnt++;
			}
		}
	}
	return (int)(cnt * SQR(resol));
  }
}

int CheckActinAbpOverlapBoundary(double *r) {
  int CS, k;
  double dist, cenDom[NDIM];

  FOR_NDIM(k) {
    cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
  }
  dist = CalcDist(cenDom, r, 0);
  CS = (dist > bnd.radRnd) ? 0 : 1;
  return CS;
}

/*---------------- Interactions between boundaries and others ----------------*/

/*-------------------------- Recording information ---------------------------*/

void RecordBoundaryForces(void) {
  
}

void RecordBoundaryLocation(void) {
  int k;
  FILE *fOut;

  // Record the locations of boundaries at the interval of "recBndLoc.prd".
  fOut = fopen(GenFileName("BndLoc"), "a");
  fprintf(fOut, "%lld\t", currTimeStep);
  FOR_NDIM(k) {
	fprintf(fOut, "%g\t%g\t", rGrid[k][0], rGrid[k][nGrid[k] - 1]);
  }
  fprintf(fOut, "\n");
  fclose(fOut);
}

void RecordBoundaryActinUnbindBind(void) {
  int n, k, *cntBndUnbRebAll, *cntBndUnbRebSum;
  FILE *fOut;

  MALLOC(cntBndUnbRebAll,int,nCpu*NDIM*4);
  MALLOC(cntBndUnbRebSum,int,NDIM*4);
  V6COPY(cntBndUnbRebSum, bndUnb.cntMe);
  V6COPY(&cntBndUnbRebSum[NDIM*2], bndReb.cntMe);

  MPI_Gather(cntBndUnbRebSum, NDIM * 4, MPI_INT, cntBndUnbRebAll, NDIM * 4, 
		MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 0; n < NDIM * 4; n++) {
		for(k = 1; k < nCpu; k++) {
			cntBndUnbRebSum[n] += P2A(cntBndUnbRebAll,k,n,NDIM*4);
		}
	}
	fOut = fopen(GenFileName("BndUnbReb"), "a");
	fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayInt(fOut, cntBndUnbRebSum, NDIM * 4, 0);
	fclose(fOut);
  }
  free(cntBndUnbRebSum);
  free(cntBndUnbRebAll);
}

void RecordBoundaryActinMature(void) {
}

void RecordBoundaryTractionForce(void) {
}


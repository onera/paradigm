/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_bound.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_part_bound_priv.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_part_bound_t *
PDM_part_bound_create
(
const int                    lComm,
const int                    nElt,
const int                    nEltPartBound,
const PDM_part_bound_cplx_t  cplx,
const int                   *nConnectedElt,
const int                   *nOfferElt,
const PDM_g_num_t            nTotalOfferElt,
const int                   nLocalOfferElt,
const PDM_g_num_t            *localOfferLnToGn
)

{
  PDM_part_bound_t * part_bound;
  _part_bound_t *tmp_part_bound;
  PDM_malloc(tmp_part_bound, 1, _part_bound_t);
  part_bound = (PDM_part_bound_t *) tmp_part_bound;
  _part_bound_t *_part_bound = (_part_bound_t *)  part_bound;

  _part_bound->nElt = nElt;
  _part_bound->nEltPartBound = nEltPartBound;
  _part_bound->lComm = lComm;
  PDM_malloc(_part_bound->eltPartBoundIdx, nEltPartBound + 1, int);
  _part_bound->eltPartBound = NULL;
  _part_bound->cplx = cplx;
  _part_bound->eltPartBoundIdx[0] = 0;
  _part_bound->nTotalOfferElt = nTotalOfferElt;
  _part_bound->nLocalOfferElt = nLocalOfferElt;
  _part_bound->localOfferLnToGn = localOfferLnToGn;

  if (cplx == PDM_PART_BOUND_SIMPLE) {
    PDM_malloc(_part_bound->nConnectedElt  , 1                , int);
    PDM_malloc(_part_bound->connectedEltIdx, nEltPartBound + 1, int);
    _part_bound->nConnectedElt[0] = *nConnectedElt;
    PDM_malloc(_part_bound->eltPartBound, (nDataEltPartBoundIni + nDataEltPartBoundElt * (*nConnectedElt)) * _part_bound->nEltPartBound, int);
    _part_bound->connectedEltIdx[0] = 0;

    PDM_malloc(_part_bound->nOfferElt  , 1                , int);
    PDM_malloc(_part_bound->offerEltIdx, nEltPartBound + 1, int);
    _part_bound->nOfferElt[0] = *nOfferElt;
    _part_bound->offerEltIdx[0] = 0;
    PDM_malloc(_part_bound->offerElt   , nEltPartBound * (*nOfferElt), int);
    PDM_malloc(_part_bound->offerLnToGn, nEltPartBound * (*nOfferElt), PDM_g_num_t);

    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->connectedEltIdx[i+1] = _part_bound->connectedEltIdx[i] +
	*nConnectedElt;
      _part_bound->offerEltIdx[i+1] =  _part_bound->offerEltIdx[i] + *nOfferElt;
      _part_bound->eltPartBoundIdx[i+1] = _part_bound->eltPartBoundIdx[i] +
                                          nDataEltPartBoundIni +
                                          nDataEltPartBoundElt * (*nConnectedElt);
    }
  }
  else {
    PDM_malloc(_part_bound->nConnectedElt  , nEltPartBound    , int);
    PDM_malloc(_part_bound->connectedEltIdx, nEltPartBound + 1, int);
    _part_bound->connectedEltIdx[0] = 0;
    memcpy(_part_bound->nConnectedElt, nConnectedElt, sizeof(int)*nEltPartBound);

    PDM_malloc(_part_bound->nOfferElt  , nEltPartBound    , int);
    PDM_malloc(_part_bound->offerEltIdx, nEltPartBound + 1, int);
    _part_bound->offerEltIdx[0] = 0;
    memcpy(_part_bound->nOfferElt, nOfferElt, sizeof(int)*nEltPartBound);

    int tConnectedElt = 0;
    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->connectedEltIdx[i+1] = _part_bound->connectedEltIdx[i] +
	nConnectedElt[i];
      _part_bound->offerEltIdx[i+1] = _part_bound->offerEltIdx[i] + nOfferElt[i];
      tConnectedElt += nConnectedElt[i];
    }

    PDM_malloc(_part_bound->offerElt   , _part_bound->offerEltIdx[nEltPartBound], int);
    PDM_malloc(_part_bound->offerLnToGn, _part_bound->offerEltIdx[nEltPartBound], PDM_g_num_t);

    PDM_malloc(_part_bound->eltPartBound, (nDataEltPartBoundIni * _part_bound->nEltPartBound + nDataEltPartBoundElt * tConnectedElt), int);

    for (int i = 0; i < nEltPartBound; i++) {
      _part_bound->eltPartBoundIdx[i+1] = _part_bound->eltPartBoundIdx[i] +
                                           nDataEltPartBoundIni +
                                           nDataEltPartBoundElt * nConnectedElt[i];
    }
  }

  PDM_malloc(_part_bound->localElt2BoundElt, nElt, int);

  return part_bound;
}

void
PDM_part_bound_local_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         localElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->eltPartBoundIdx[iBoundElt];
  _part_bound->eltPartBound[idx] = localElt;
  _part_bound->localElt2BoundElt[localElt - 1] = boundElt;
}


int
PDM_part_bound_n_offer_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int nOffer;
  if (_part_bound->cplx == PDM_PART_BOUND_CPLX) {
    nOffer = _part_bound->nOfferElt[iBoundElt];
  }
  else {
    nOffer = _part_bound->nOfferElt[0];
  }
  return nOffer;
}


PDM_g_num_t
PDM_part_bound_n_total_offer_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nTotalOfferElt;
}


int
PDM_part_bound_n_local_offer_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nLocalOfferElt;
}


const PDM_g_num_t *
PDM_part_bound_local_offer_elt_ln_to_gn_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->localOfferLnToGn;
}


void
PDM_part_bound_offer_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iOfferElt,
 const int         lNum,
 const PDM_g_num_t  gNum
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->offerEltIdx[iBoundElt] + iOfferElt;
  _part_bound->offerElt[idx] = lNum;
  _part_bound->offerLnToGn[idx] = gNum;
}

void
PDM_part_bound_offer_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iOfferElt,
 int              *lNum,
 PDM_g_num_t       *gNum
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->offerEltIdx[iBoundElt] + iOfferElt;
  *lNum = _part_bound->offerElt[idx];
  *gNum = _part_bound->offerLnToGn[idx];
}


int
PDM_part_bound_n_elt_get
(
 PDM_part_bound_t *part_bound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  return _part_bound->nElt;
}

int
PDM_part_bound_n_elt_bound_get
(
 PDM_part_bound_t *part_bound
)
{
 _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 return _part_bound->nEltPartBound;
}


PDM_part_bound_cplx_t
PDM_part_bound_cplx_get
(
 PDM_part_bound_t *part_bound
)
{
 _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 return _part_bound->cplx;
}


void
PDM_part_bound_bound_elt_get
(
 PDM_part_bound_t *part_bound,
 const int      boundElt,
       int     *localElt,
       int     *nConnectedElt
 )
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;
  int idx = _part_bound->eltPartBoundIdx[iBoundElt];
  *localElt = _part_bound->eltPartBound[idx];
  if (_part_bound->cplx == PDM_PART_BOUND_CPLX) {
    *nConnectedElt = _part_bound->nConnectedElt[iBoundElt];
  }
  else {
    *nConnectedElt = _part_bound->nConnectedElt[0];
  }
}


void
PDM_part_bound_local_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         localElt,
       int        *boundElt,
       int        *nConnectedElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iLocalElt = localElt - 1;
  *boundElt = _part_bound->localElt2BoundElt[iLocalElt];

  int localElt2;
  PDM_part_bound_bound_elt_get (part_bound,
                                *boundElt,
                                &localElt2,
                                nConnectedElt);

}


void
PDM_part_bound_distant_elt_set
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iConnectedElt,
 const int         iProc,
 const int         iProcPart,
 const int         iProcPartElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;

  int nConnectedElt = _part_bound->nConnectedElt[0];

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    nConnectedElt = _part_bound->nConnectedElt[iBoundElt];

  if (iConnectedElt >= nConnectedElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error part_bound_distant_elt_set :"
	    "Error in edgeFace computing\n");
    abort();
  }

  int idx = _part_bound->eltPartBoundIdx[iBoundElt] + nDataEltPartBoundIni +
    iConnectedElt * nDataEltPartBoundElt;

  _part_bound->eltPartBound[idx++] = iProc;
  _part_bound->eltPartBound[idx++] = iProcPart;
  _part_bound->eltPartBound[idx++] = iProcPartElt;
  _part_bound->eltPartBound[idx++] = _part_bound->connectedEltIdx[iBoundElt] +
    iConnectedElt;

}



void
PDM_part_bound_distant_elt_get
(
 PDM_part_bound_t *part_bound,
 const int         boundElt,
 const int         iConnectedElt,
       int        *iProc,
       int        *iProcPart,
       int        *iProcPartElt,
       int        *iDistElt
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  int iBoundElt = boundElt - 1;

  int nConnectedElt = _part_bound->nConnectedElt[0];

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    nConnectedElt = _part_bound->nConnectedElt[iBoundElt];

  if (iConnectedElt >= nConnectedElt) {
    PDM_error(__FILE__, __LINE__, 0, "Error part_bound_distant_elt_get :"
	    "iConnectedElt > nConnectedElt\n");
    abort();
  }

  int idx = _part_bound->eltPartBoundIdx[iBoundElt] +
    nDataEltPartBoundIni +
    iConnectedElt * nDataEltPartBoundElt;

  *iProc        = _part_bound->eltPartBound[idx++];
  *iProcPart    = _part_bound->eltPartBound[idx++];
  *iProcPartElt = _part_bound->eltPartBound[idx++];
  *iDistElt     = _part_bound->eltPartBound[idx++];
}


void
PDM_part_bound_adjust_size
(
 PDM_part_bound_t *part_bound,
 const int         nEltPartBound
)
{
  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
  if (nEltPartBound > _part_bound->nEltPartBound) {
    PDM_error(__FILE__, __LINE__, 0, "Error _part_bound_adjust_size : Error this function"
                    "can't increase the size of part_bound structure\n");
    abort();
  }

  _part_bound->nEltPartBound = nEltPartBound;
  PDM_realloc(_part_bound->eltPartBoundIdx ,_part_bound->eltPartBoundIdx , (nEltPartBound + 1),int);

  if (_part_bound->cplx == PDM_PART_BOUND_CPLX)
    PDM_realloc(_part_bound->nConnectedElt ,_part_bound->nConnectedElt , nEltPartBound,int);

  PDM_realloc(_part_bound->eltPartBound, _part_bound->eltPartBound, _part_bound->eltPartBoundIdx[nEltPartBound], int);

}



PDM_part_bound_t *
PDM_part_bound_free
(
PDM_part_bound_t *part_bound
)
{

  _part_bound_t *_part_bound = (_part_bound_t *) part_bound;
 if (_part_bound != NULL) {
    if (_part_bound->eltPartBoundIdx != NULL)
      PDM_free(_part_bound->eltPartBoundIdx);
    if (_part_bound->eltPartBound != NULL)
      PDM_free(_part_bound->eltPartBound);
    if (_part_bound->connectedEltIdx != NULL)
      PDM_free(_part_bound->connectedEltIdx);
    if (_part_bound->nConnectedElt != NULL)
      PDM_free(_part_bound->nConnectedElt);
    if (_part_bound->localElt2BoundElt != NULL)
      PDM_free(_part_bound->localElt2BoundElt);
    if (_part_bound->nOfferElt != NULL)
      PDM_free(_part_bound->nOfferElt);
    if (_part_bound->offerEltIdx != NULL)
      PDM_free(_part_bound->offerEltIdx);
    if (_part_bound->offerElt != NULL)
      PDM_free(_part_bound->offerElt);
    if (_part_bound->offerLnToGn != NULL)
      PDM_free(_part_bound->offerLnToGn);

    PDM_free(_part_bound);
  }
  return NULL;
}



void
PDM_part_bound_dump
(
 PDM_part_bound_t *part_bound
 )
{
  if (part_bound != NULL) {

    PDM_part_bound_cplx_t cplx = PDM_part_bound_cplx_get (part_bound);
    int nEltPartBound = PDM_part_bound_n_elt_bound_get (part_bound);
    PDM_printf("PDM_part_bound :\n");
    PDM_printf("  - cplx %d\n", cplx);
    PDM_printf("  - nEltPartBound %d\n", nEltPartBound);
    for (int j = 0; j < nEltPartBound; j++) {

      int nConnectedElt;
      int localElt;

      PDM_part_bound_bound_elt_get (part_bound,
                                    j+1,
                                    &localElt,
                                    &nConnectedElt);

      PDM_printf("    - localElt %d\n", localElt);
      PDM_printf("      - nConnectedElt %d\n", nConnectedElt);
      for (int k = 0; k < nConnectedElt; k++) {
        int iProc;
        int i_part;
        int iElt;
        int iDistElt;
        PDM_part_bound_distant_elt_get (part_bound,
                                        j+1,
                                        k,
                                        &iProc,
                                        &i_part,
                                        &iElt,
                                        &iDistElt);
        PDM_printf("        - %d %d %d %d (iproc i_part, iElt, iDistElt)\n", iProc, i_part, iElt, iDistElt);

      }
      int nOfferElt = PDM_part_bound_n_offer_elt_get (part_bound, j+1);
      PDM_printf("      - nOfferElt %d\n", nOfferElt);

      for (int k = 0; k < nOfferElt; k++) {
        int        lNum;
        PDM_g_num_t gNum;
        PDM_part_bound_offer_elt_get (part_bound,
                                      j+1,
                                      k,
                                      &lNum,
                                      &gNum);
        PDM_printf("        - %d "PDM_FMT_G_NUM" (lNum gNum)\n", lNum, gNum);

      }
    }

    fflush(stdout);
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

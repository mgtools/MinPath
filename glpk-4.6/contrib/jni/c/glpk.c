/* glpk.c (Java Native Interface for GLPK) */

/* ---------------------------------------------------------------------
-- Copyright (C) 2003 Andrew Makhorin <mao@mai2.rcnet.ru>, Department
-- for Applied Informatics, Moscow Aviation Institute, Moscow, Russia.
-- All rights reserved.
--
-- Author: Yuri Victorovich, Software Engineer, yuri@gjt.org.
--  
-- This file is a part of GLPK (GNU Linear Programming Kit).
--  
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--  
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--  
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
-- 02111-1307, USA.
------------------------------------------------------------------------*/


#include <jni.h>
#include <glpk.h>
#include <malloc.h>

#include "org_gnu_glpk_GlpkSolver.h"

#define CLS_NAME "org/gnu/glpk/GlpkSolver"

static void put_lpx(JNIEnv *env, jobject obj, LPX *lp) {
  jclass cls = (*env)->GetObjectClass(env, obj);
  jfieldID fid = (*env)->GetFieldID(env, cls, "lp", "I");
  (*env)->SetIntField(env, obj, fid, (jint)lp);
}

static LPX * get_lpx(JNIEnv *env, jobject obj) {
  LPX *lp;

  jclass cls = (*env)->GetObjectClass(env, obj);
  jfieldID fid = (*env)->GetFieldID(env, cls, "lp", "I");
  lp = (LPX*)(*env)->GetIntField(env, obj, fid);
  if (lp != 0)
    return (lp);

  lp = lpx_create_prob();
  put_lpx(env, obj, lp);
  return (lp);
}

static jobject create_string(JNIEnv *env, char *str) {
  return (*env)->NewStringUTF(env, str);
}


JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_finalize(JNIEnv *env, jobject obj) {
  LPX *lp = get_lpx(env, obj);
  if (lp != 0) {
    lpx_delete_prob(lp);
    put_lpx(env, obj, 0);
  }
}

typedef struct {
    JNIEnv  *env;
    jobject obj;
} info_t;

int _hook_fn(info_t *info, char *msg, int bFault) {
	int debug=0;

    jobject obj=info->obj;
    JNIEnv  *env=info->env;

	if (debug) printf("HERE0: %s\n",msg);

	jclass cls=(*env)->FindClass(env,"org/gnu/glpk/GlpkSolver");

    if (cls == NULL)
		return 0;

	if (debug) printf("HERE1\n");

	jstring	callbackArg=(*env)->NewStringUTF(env, msg);

	if (debug) printf("HERE2\n");

    jmethodID 	mID;

    if (bFault)
        mID = (*env)->GetMethodID(env,cls,"faultHook","(Ljava/lang/String;)V");
	else
        mID = (*env)->GetMethodID(env,cls,"printHook","(Ljava/lang/String;)V");

    if (mID == NULL)
        return 0;

    if (debug) printf("HERE3\n");

    (*env)->CallVoidMethod(env, obj, mID, callbackArg);

	if (debug) printf("HERE4\n");

	return 1;
}

int _fault_hook_fn(void *info, char *msg) {
	return _hook_fn((info_t*)info,msg,1);
}

int _print_hook_fn(void *info, char *msg) {
    return _hook_fn((info_t*)info,msg,0);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_enablePrints(JNIEnv *env, jobject obj, jboolean enable)
{
    info_t  *info=malloc(sizeof(info_t));

    info->env = env;
    info->obj = obj;

    lib_set_fault_hook((void*)info, enable ? NULL : &_fault_hook_fn);
    lib_set_print_hook((void*)info, enable ? NULL : &_print_hook_fn);
}

/* Problem creating and modifying routines */

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_addRows(JNIEnv *env, jobject obj, jint nrs) {
  lpx_add_rows(get_lpx(env, obj), (int)nrs);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_addCols(JNIEnv *env, jobject obj, jint ncs) {
  lpx_add_cols(get_lpx(env, obj), (int)ncs);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_checkName(JNIEnv *env, jclass cls, jstring name) {
  const char* chr = (*env)->GetStringUTFChars(env, name, 0);
  jboolean ret = lpx_check_name((char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, name, chr);

  return  (ret);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setProbName(JNIEnv *env, jobject obj,
  jstring name)
{
  LPX *lp = get_lpx(env, obj);
  char * chr = (char *)(*env)->GetStringUTFChars(env, name, 0);
  lpx_set_prob_name(lp, chr);
  (*env)->ReleaseStringUTFChars(env, name, chr);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setRowName(JNIEnv *env, jobject obj,
  jint i, jstring name)
{
  LPX *lp = get_lpx(env, obj);
  lpx_set_row_name(lp, i, (char*)(*env)->GetStringUTFChars(env, name, 0));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setColName(JNIEnv *env, jobject obj,
  jint j, jstring name)
{
  LPX *lp = get_lpx(env, obj);
  lpx_set_col_name(lp, j, (char*)(*env)->GetStringUTFChars(env, name, 0));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setRowBnds(JNIEnv *env, jobject obj, jint i,
  jint typx, jdouble lb, jdouble ub)
{
  LPX *lp = get_lpx(env, obj);
  lpx_set_row_bnds(lp, (int)i, (int)typx, (double)lb, (double)ub);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setColBnds(JNIEnv *env, jobject obj, jint j,
  jint typx, jdouble lb, jdouble ub)
{
  LPX *lp = get_lpx(env, obj);
  lpx_set_col_bnds(lp, (int)j, (int)typx, (double)lb, (double)ub);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setObjName(JNIEnv *env, jobject obj,
  jstring name)
{
  LPX *lp = get_lpx(env, obj);
  lpx_set_obj_name(lp, (char*)(*env)->GetStringUTFChars(env, name, 0));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setObjDir(JNIEnv *env, jobject obj, jint dir) {
  lpx_set_obj_dir(get_lpx(env, obj), (int)dir);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setObjC0(JNIEnv *env, jobject obj, jdouble c0) {
  lpx_set_obj_dir(get_lpx(env, obj), (double)c0);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setRowCoef(JNIEnv *env, jobject obj,
  jint i, jdouble coef)
{
  lpx_set_row_coef(get_lpx(env, obj), (int)i, (double)coef);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setColCoef(JNIEnv *env, jobject obj,
  jint j, jdouble coef)
{
  lpx_set_col_coef(get_lpx(env, obj), (int)j, (double)coef);
}

struct env_cls {
  JNIEnv *env;
  jclass loaderObj;
  jmethodID matMid;
  jfieldID iFid;
  jfieldID jFid;
  jfieldID coefFid;
  jobject elem;
};

static double local_mat_fn(void *info, int *i, int *j) {
  double ret;
  struct env_cls *ec = (struct env_cls *)info;

  if (!(*ec->env)->CallBooleanMethod(ec->env, ec->loaderObj, ec->matMid,
      ec->elem)) { // false returned: end of file
    *i = *j = 0;
    return (0.0);
  }
  if ((*ec->env)->ExceptionOccurred(ec->env)) { // some exception occured
    *i = *j = 0;
    return (0);
  }

  *i = (int)(*ec->env)->GetIntField(ec->env, ec->elem, ec->iFid);
  *j = (int)(*ec->env)->GetIntField(ec->env, ec->elem, ec->jFid);
  ret = (double)(*ec->env)->GetDoubleField(ec->env, ec->elem, ec->coefFid);

  return (ret);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_loadMat(JNIEnv *env, jobject obj, jobject loader) {
  LPX *lp = get_lpx(env, obj);
  jclass loaderCls, elementCls;
  jmethodID elementMid;
  struct env_cls ec;

  loaderCls = (*env)->GetObjectClass(env, loader);
  elementCls = (*env)->FindClass(env, CLS_NAME "MatElement");
  if (elementCls == 0)
    return;

  ec.env = env;
  ec.loaderObj = loader;
  if ((ec.matMid = (*env)->GetMethodID(env, loaderCls, "mat",
                      "(L" CLS_NAME "MatElement;)Z")) == 0 ||
      (elementMid = (*env)->GetMethodID(env, elementCls, "<init>",
                      "()V")) == 0 ||
      (ec.iFid    = (*env)->GetFieldID(env, elementCls, "i", "I")) == 0 ||
      (ec.jFid    = (*env)->GetFieldID(env, elementCls, "j", "I")) == 0 ||
      (ec.coefFid = (*env)->GetFieldID(env, elementCls, "coef", "D")) == 0)
    return;

  ec.elem = (*env)->NewObject(env, elementCls, elementMid);
  if (ec.elem == 0)
    return; /* exception occured */

  lpx_load_mat(lp, &ec, local_mat_fn);
  // if exception occured it should be carried through DeleteLocalRef

  (*env)->DeleteLocalRef(env, ec.elem);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_loadMat3(JNIEnv *env, jobject obj,
  jint nz, jintArray rn, jintArray cn, jdoubleArray a)
{
  LPX *lp = get_lpx(env, obj);
  jint *rnEls   = (*env)->GetIntArrayElements(env, rn, 0);
  jint *cnEls   = (*env)->GetIntArrayElements(env, cn, 0);
  jdouble *aEls = (*env)->GetDoubleArrayElements(env, a, 0);

  lpx_load_mat3(lp, (int)nz, (int*)rnEls, (int*)cnEls, (double*)aEls);

  (*env)->ReleaseIntArrayElements(env, rn, rnEls, 0);
  (*env)->ReleaseIntArrayElements(env, cn, cnEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, cn, aEls, 0);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setMatRow(JNIEnv *env, jobject obj,
  jint i, jint len, jintArray ndx, jdoubleArray val)
{
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  lpx_set_mat_row(lp, (int)i, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);
}


JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setMatCol(JNIEnv *env, jobject obj,
  jint j, jint len, jintArray ndx, jdoubleArray val)
{
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  lpx_set_mat_col(lp, (int)j, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_unmarkAll(JNIEnv *env, jobject obj) {
  LPX *lp = get_lpx(env, obj);
  lpx_unmark_all(lp);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_markRow(JNIEnv *env, jobject obj, jint i, jint mark) {
  LPX *lp = get_lpx(env, obj);
  lpx_mark_row(lp, (int)i, (int)mark);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_markCol(JNIEnv *env, jobject obj, jint j, jint mark) {
  LPX *lp = get_lpx(env, obj);
  lpx_mark_col(lp, (int)j, (int)mark);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_clearMat(JNIEnv *env, jobject obj) {
  lpx_clear_mat(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_delItems(JNIEnv *env, jobject obj) {
  lpx_del_items(get_lpx(env, obj));
}

/* Problem querying routines */

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getNumRows(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_num_rows(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getNumCols(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_num_cols(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getNumNz(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_num_nz(get_lpx(env, obj));
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_getProbName(JNIEnv *env, jobject obj) {
  return create_string(env, lpx_get_prob_name(get_lpx(env, obj)));
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_getRowName(JNIEnv *env, jobject obj, jint i) {
  return create_string(env, lpx_get_row_name(get_lpx(env, obj), (int)i));
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_getColName(JNIEnv *env, jobject obj, jint j) {
  return create_string(env, lpx_get_col_name(get_lpx(env, obj), (int)j));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getRowBnds(JNIEnv *env, jobject obj,
  jint i, jobject bnds)
{
  LPX *lp = get_lpx(env, obj);
  jclass cls = (*env)->GetObjectClass(env, bnds);
  int typx;
  double lb, ub;
  jfieldID fidTypx, fidLB, fidUB;
  if ((fidTypx = (*env)->GetFieldID(env, cls, "typx", "I")) == 0 ||
      (fidLB   = (*env)->GetFieldID(env, cls, "lb", "D")) == 0 ||
      (fidUB   = (*env)->GetFieldID(env, cls, "ub", "D")) == 0)
    return;

  lpx_get_row_bnds(lp, (int)i, &typx, &lb, &ub);

  (*env)->SetIntField(env, bnds, fidTypx, (jint)typx);
  (*env)->SetDoubleField(env, bnds, fidLB, (jdouble)lb);
  (*env)->SetDoubleField(env, bnds, fidUB, (jdouble)ub);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getColBnds(JNIEnv *env, jobject obj, jint j,
  jobject bnds)
{
  LPX *lp = get_lpx(env, obj);
  jclass cls = (*env)->GetObjectClass(env, bnds);
  int typx;
  double lb, ub;
  jfieldID fidTypx, fidLB, fidUB;
  if ((fidTypx = (*env)->GetFieldID(env, cls, "typx", "I")) == 0 ||
      (fidLB   = (*env)->GetFieldID(env, cls, "lb", "D")) == 0 ||
      (fidUB   = (*env)->GetFieldID(env, cls, "ub", "D")) == 0)
    return;

  lpx_get_col_bnds(lp, (int)j, &typx, &lb, &ub);

  (*env)->SetIntField(env, bnds, fidTypx, (jint)typx);
  (*env)->SetDoubleField(env, bnds, fidLB, (jdouble)lb);
  (*env)->SetDoubleField(env, bnds, fidUB, (jdouble)ub);
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_getObjName(JNIEnv *env, jobject obj) {
  return create_string(env, lpx_get_obj_name(get_lpx(env, obj)));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getObjDir(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_obj_dir(get_lpx(env, obj));
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getObjC0(JNIEnv *env, jobject obj) {
  return (jdouble)lpx_get_obj_c0(get_lpx(env, obj));
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getRowCoef(JNIEnv *env, jobject obj, jint i) {
  return (jdouble)lpx_get_row_coef(get_lpx(env, obj), (int)i);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getColCoef(JNIEnv *env, jobject obj, jint j) {
  return (jdouble)lpx_get_col_coef(get_lpx(env, obj), (int)j);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getMatRow(JNIEnv *env, jobject obj, jint i,
  jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_get_mat_row(lp, (int)i, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return ((jint)ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getMatCol(JNIEnv *env, jobject obj, jint j,
  jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_get_mat_col(lp, (int)j, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return ((jint)ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getRowMark(JNIEnv *env, jobject obj, jint i) {
  return (jint)lpx_get_row_mark(get_lpx(env, obj), (int)i);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getColMark(JNIEnv *env, jobject obj, jint j) {
  return (jint)lpx_get_col_mark(get_lpx(env, obj), (int)j);
}

/* Problem scaling routines */

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_scaleProb(JNIEnv *env, jobject obj) {
  lpx_scale_prob(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_unscaleProb(JNIEnv *env, jobject obj) {
  lpx_unscale_prob(get_lpx(env, obj));
}

/* Basis constructing routines */

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_stdBasis(JNIEnv *env, jobject obj) {
  lpx_std_basis(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_advBasis(JNIEnv *env, jobject obj) {
  lpx_adv_basis(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setRowStat(JNIEnv *env, jobject obj,
  jint i, jint stat)
{
  lpx_set_row_stat(get_lpx(env, obj), (int)i, (int)stat);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setColStat(JNIEnv *env, jobject obj,
  jint j, jint stat)
{
  lpx_set_col_stat(get_lpx(env, obj), (int)j, (int)stat);
}

/* Simplex method routines */

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_warmUp(JNIEnv *env, jobject obj) {
  return (jint)lpx_warm_up(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_simplex(JNIEnv *env, jobject obj) {
  return (jint)lpx_simplex(get_lpx(env, obj));
}

/* Basic solutin querying routines */

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getStatus(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_status(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getPrimStat(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_prim_stat(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getDualStat(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_dual_stat(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getRowInfo(JNIEnv *env, jobject obj, jint i,
  jobject info)
{
  jclass cls = (*env)->GetObjectClass(env, info);
  jfieldID fidTagx, fidVx, fidDx;
  LPX *lp = get_lpx(env, obj);
  int tagx;
  double vx, dx;

  if ((fidTagx = (*env)->GetFieldID(env, cls, "tagx", "I")) == 0 ||
      (fidVx   = (*env)->GetFieldID(env, cls, "vx", "D")) == 0 ||
      (fidDx   = (*env)->GetFieldID(env, cls, "dx", "D")) == 0)
    return;

  lpx_get_row_info(lp, (int)i, &tagx, &vx, &dx);

  (*env)->SetIntField   (env, info, fidTagx, (jint)tagx);
  (*env)->SetDoubleField(env, info, fidVx,   (jdouble)vx);
  (*env)->SetDoubleField(env, info, fidDx,   (jdouble)dx);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getColInfo(JNIEnv *env, jobject obj, jint j,
  jobject info)
{
  jclass cls = (*env)->GetObjectClass(env, info);
  jfieldID fidTagx, fidVx, fidDx;
  LPX *lp = get_lpx(env, obj);
  int tagx;
  double vx, dx;

  if ((fidTagx = (*env)->GetFieldID(env, cls, "tagx", "I")) == 0 ||
      (fidVx   = (*env)->GetFieldID(env, cls, "vx", "D")) == 0 ||
      (fidDx   = (*env)->GetFieldID(env, cls, "dx", "D")) == 0)
    return;

  lpx_get_col_info(lp, (int)j, &tagx, &vx, &dx);

  (*env)->SetIntField   (env, info, fidTagx, (jint)tagx);
  (*env)->SetDoubleField(env, info, fidVx,   (jdouble)vx);
  (*env)->SetDoubleField(env, info, fidDx,   (jdouble)dx);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getObjVal(JNIEnv *env, jobject obj) {
  return (jdouble)lpx_get_obj_val(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_checkKkt(JNIEnv *env, jobject obj, jint scaled,
  jobject jkkt)
{
  LPX *lp = get_lpx(env, obj);
  LPXKKT kkt;
  jclass cls = (*env)->GetObjectClass(env, jkkt);
  jfieldID fid_pe_ae_max, fid_pe_ae_row, fid_pe_re_max, fid_pe_re_row,
           fid_pe_quality, fid_pb_ae_max, fid_pb_ae_ind, fid_pb_re_max,
           fid_pb_re_ind, fid_pb_quality, fid_de_ae_max, fid_de_ae_col,
           fid_de_re_max, fid_de_re_col, fid_de_quality, fid_db_ae_max,
           fid_db_ae_ind, fid_db_re_max, fid_db_re_ind, fid_db_quality,
           fid_cs_ae_max, fid_cs_ae_ind, fid_cs_re_max, fid_cs_re_ind,
           fid_cs_quality;
  if ((fid_pe_ae_max = (*env)->GetFieldID(env, cls, "pe_ae_max", "D")) == 0 ||
      (fid_pe_ae_row = (*env)->GetFieldID(env, cls, "pe_ae_row", "I")) == 0 ||
      (fid_pe_re_max = (*env)->GetFieldID(env, cls, "pe_re_max", "D")) == 0 ||
      (fid_pe_re_row = (*env)->GetFieldID(env, cls, "pe_re_row", "I")) == 0 ||
      (fid_pe_quality= (*env)->GetFieldID(env, cls, "pe_quality","I")) == 0 ||
      (fid_pb_ae_max = (*env)->GetFieldID(env, cls, "pb_ae_max", "D")) == 0 ||
      (fid_pb_ae_ind = (*env)->GetFieldID(env, cls, "pb_ae_ind", "I")) == 0 ||
      (fid_pb_re_max = (*env)->GetFieldID(env, cls, "pb_re_max", "D")) == 0 ||
      (fid_pb_re_ind = (*env)->GetFieldID(env, cls, "pb_re_ind", "I")) == 0 ||
      (fid_pb_quality= (*env)->GetFieldID(env, cls, "pb_quality","I")) == 0 ||
      (fid_de_ae_max = (*env)->GetFieldID(env, cls, "de_ae_max", "D")) == 0 ||
      (fid_de_ae_col = (*env)->GetFieldID(env, cls, "de_ae_col", "I")) == 0 ||
      (fid_de_re_max = (*env)->GetFieldID(env, cls, "de_re_max", "D")) == 0 ||
      (fid_de_re_col = (*env)->GetFieldID(env, cls, "de_re_col", "I")) == 0 ||
      (fid_de_quality= (*env)->GetFieldID(env, cls, "de_quality","I")) == 0 ||
      (fid_db_ae_max = (*env)->GetFieldID(env, cls, "db_ae_max", "D")) == 0 ||
      (fid_db_ae_ind = (*env)->GetFieldID(env, cls, "db_ae_ind", "I")) == 0 ||
      (fid_db_re_max = (*env)->GetFieldID(env, cls, "db_re_max", "D")) == 0 ||
      (fid_db_re_ind = (*env)->GetFieldID(env, cls, "db_re_ind", "I")) == 0 ||
      (fid_db_quality= (*env)->GetFieldID(env, cls, "db_quality","I")) == 0 ||
      (fid_cs_ae_max = (*env)->GetFieldID(env, cls, "cs_ae_max", "D")) == 0 ||
      (fid_cs_ae_ind = (*env)->GetFieldID(env, cls, "cs_ae_ind", "I")) == 0 ||
      (fid_cs_re_max = (*env)->GetFieldID(env, cls, "cs_re_max", "D")) == 0 ||
      (fid_cs_re_ind = (*env)->GetFieldID(env, cls, "cs_re_ind", "I")) == 0 ||
      (fid_cs_quality= (*env)->GetFieldID(env, cls, "cs_quality","I")) == 0)
    return;

  lpx_check_kkt(lp, (int)scaled, &kkt);

  (*env)->SetDoubleField(env, jkkt, fid_pe_ae_max,  (jdouble)kkt.pe_ae_max);
  (*env)->SetIntField   (env, jkkt, fid_pe_ae_row,  (jint)   kkt.pe_ae_row);
  (*env)->SetDoubleField(env, jkkt, fid_pe_re_max,  (jdouble)kkt.pe_re_max);
  (*env)->SetIntField   (env, jkkt, fid_pe_re_row,  (jint)   kkt.pe_re_row);
  (*env)->SetIntField   (env, jkkt, fid_pe_quality, (jint)   kkt.pe_quality);
  (*env)->SetDoubleField(env, jkkt, fid_pb_ae_max,  (jdouble)kkt.pb_ae_max);
  (*env)->SetIntField   (env, jkkt, fid_pb_ae_ind,  (jint)   kkt.pb_ae_ind);
  (*env)->SetDoubleField(env, jkkt, fid_pb_re_max,  (jdouble)kkt.pb_re_max);
  (*env)->SetIntField   (env, jkkt, fid_pb_re_ind,  (jint)   kkt.pb_re_ind);
  (*env)->SetIntField   (env, jkkt, fid_pb_quality, (jint)   kkt.pb_quality);
  (*env)->SetDoubleField(env, jkkt, fid_de_ae_max,  (jdouble)kkt.de_ae_max);
  (*env)->SetIntField   (env, jkkt, fid_de_ae_col,  (jint)   kkt.de_ae_col);
  (*env)->SetDoubleField(env, jkkt, fid_de_re_max,  (jdouble)kkt.de_re_max);
  (*env)->SetIntField   (env, jkkt, fid_de_re_col,  (jint)   kkt.de_re_col);
  (*env)->SetIntField   (env, jkkt, fid_de_quality, (jint)   kkt.de_quality);
  (*env)->SetDoubleField(env, jkkt, fid_db_ae_max,  (jdouble)kkt.db_ae_max);
  (*env)->SetIntField   (env, jkkt, fid_db_ae_ind,  (jint)   kkt.db_ae_ind);
  (*env)->SetDoubleField(env, jkkt, fid_db_re_max,  (jdouble)kkt.db_re_max);
  (*env)->SetIntField   (env, jkkt, fid_db_re_ind,  (jint)   kkt.db_re_ind);
  (*env)->SetIntField   (env, jkkt, fid_db_quality, (jint)   kkt.db_quality);
  (*env)->SetDoubleField(env, jkkt, fid_cs_ae_max,  (jdouble)kkt.cs_ae_max);
  (*env)->SetIntField   (env, jkkt, fid_cs_ae_ind,  (jint)   kkt.cs_ae_ind);
  (*env)->SetDoubleField(env, jkkt, fid_cs_re_max,  (jdouble)kkt.cs_re_max);
  (*env)->SetIntField   (env, jkkt, fid_cs_re_ind,  (jint)   kkt.cs_re_ind);
  (*env)->SetIntField   (env, jkkt, fid_cs_quality, (jint)   kkt.cs_quality);
}


/* Simplex table routines */

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_evalTabRow(JNIEnv *env, jobject obj,
  jint k, jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_eval_tab_row(lp, (int)k, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_evalTabCol(JNIEnv *env, jobject obj,
  jint k, jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_eval_tab_col(lp, (int)k, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_transformRow(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_transform_row(lp, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_transformCol(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_transform_col(lp, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_primRatioTest(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val, jint how, jdouble tol)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_prim_ratio_test(lp, (int)len, (int*)ndxEls, (double*)valEls,
                                (int)how, (double)tol);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_dualRatioTest(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val, jint how, jdouble tol)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_dual_ratio_test(lp, (int)len, (int*)ndxEls, (double*)valEls,
                                (int)how, (double)tol);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return ((jint)ret);
}

/* Interior point method routines */

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_interior(JNIEnv *env, jobject obj) {
  return (jint)lpx_interior(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getIpsStat(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_ips_stat(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getIpsRow(JNIEnv *env, jobject obj,
  jint i, jobject point)
{
  LPX *lp = get_lpx(env, obj);
  double vx, dx;
  jclass cls = (*env)->GetObjectClass(env, point);
  jfieldID fidVx, fidDx;
  if ((fidVx = (*env)->GetFieldID(env, cls, "vx", "D")) == 0 ||
      (fidDx = (*env)->GetFieldID(env, cls, "dx", "D")) == 0)
    return;

  lpx_get_ips_row(lp, (int)i, &vx, &dx);

  (*env)->SetDoubleField(env, point, fidVx, (jdouble)vx);
  (*env)->SetDoubleField(env, point, fidDx, (jdouble)dx);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_getIpsCol(JNIEnv *env, jobject obj,
  jint j, jobject point)
{
  LPX *lp = get_lpx(env, obj);
  double vx, dx;
  jclass cls = (*env)->GetObjectClass(env, point);
  jfieldID fidVx, fidDx;
  if ((fidVx = (*env)->GetFieldID(env, cls, "vx", "D")) == 0 ||
      (fidDx = (*env)->GetFieldID(env, cls, "dx", "D")) == 0)
    return;

  lpx_get_ips_col(lp, (int)j, &vx, &dx);

  (*env)->SetDoubleField(env, point, fidVx, (jdouble)vx);
  (*env)->SetDoubleField(env, point, fidDx, (jdouble)dx);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getIpsObj(JNIEnv *env, jobject obj) {
  return (jdouble)lpx_get_ips_obj(get_lpx(env, obj));
}

/* MIP routines */
JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setClss(JNIEnv *env, jobject obj, jint clss) {
  lpx_set_class(get_lpx(env, obj), (int)clss);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getClss(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_class(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setColKind(JNIEnv *env, jobject obj,
  jint j, jint kind)
{
  lpx_set_col_kind(get_lpx(env, obj), (int)j, (int)kind);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getColKind(JNIEnv *env, jobject obj, jint j) {
  return (jint)lpx_get_col_kind(get_lpx(env, obj), (int)j);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getNumInt(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_num_int(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getNumBin(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_num_bin(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_integer(JNIEnv *env, jobject obj) {
  return (jint)lpx_integer(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getMIPStat(JNIEnv *env, jobject obj) {
  return (jint)lpx_get_mip_stat(get_lpx(env, obj));
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getMIPRow(JNIEnv *env, jobject obj, jint i) {
  return (jdouble)lpx_get_mip_row(get_lpx(env, obj), (int)i);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getMIPCol(JNIEnv *env, jobject obj, jint j) {
  return (jdouble)lpx_get_mip_col(get_lpx(env, obj), (int)j);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getMIPObj(JNIEnv *env, jobject obj) {
  return (jdouble)lpx_get_mip_obj(get_lpx(env, obj));
}

/* Control parameters and statistics routines */

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_resetParms(JNIEnv *env, jobject obj) {
  lpx_reset_parms(get_lpx(env, obj));
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setIntParm(JNIEnv *env, jobject obj,
  jint param, jint val)
{
  lpx_set_int_parm(get_lpx(env, obj), (int)param, (int)val);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_getIntParm(JNIEnv *env, jobject obj, jint param) {
  return (jint)lpx_get_int_parm(get_lpx(env, obj), (int)param);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_setRealParm(JNIEnv *env, jobject obj,
  jint param, jdouble val)
{
  lpx_set_real_parm(get_lpx(env, obj), (int)param, (double)val);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_getRealParm(JNIEnv *env, jobject obj, jint param) {
  return (jdouble)lpx_get_real_parm(get_lpx(env, obj), (int)param);
}

/* Utility routines */

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_readMps(JNIEnv *env, jclass cls, jstring fname) {
  jobject ret;
  jmethodID mid;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  LPX *lp = lpx_read_mps((char*)chr);
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  if (lp == 0)
    return (0);

  mid = (*env)->GetMethodID(env, cls, "<init>", "()V");
  if (mid == 0)  {
    lpx_delete_prob(lp);
    return (0);
  }
  ret = (*env)->NewObject(env, cls, mid);
  if ((*env)->ExceptionOccurred(env)) {
    lpx_delete_prob(lp);
    return (0);
  }

  put_lpx(env, ret, lp);
  return (ret);
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_readLpt(JNIEnv *env, jclass cls, jstring fname) {
  jobject ret;
  jmethodID mid;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  LPX *lp = lpx_read_lpt((char*)chr);
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  if (lp == 0)
    return (0);

  mid = (*env)->GetMethodID(env, cls, "<init>", "()V");
  if (mid == 0) {
    lpx_delete_prob(lp);
    return (0);
  }
  ret = (*env)->NewObject(env, cls, mid);
  if ((*env)->ExceptionOccurred(env)) {
    lpx_delete_prob(lp);
    return (0);
  }

  put_lpx(env, ret, lp);
  return (ret);
}

JNIEXPORT jobject JNICALL
Java_org_gnu_glpk_GlpkSolver_readModel(JNIEnv *env, jclass cls,
  jstring model, jstring data, jstring output)
{
  jobject ret;
  jmethodID mid;
  LPX *lp;
  const char *chr_model = 0, *chr_data, *chr_output  = 0;

  if (model != 0)
    chr_model = (*env)->GetStringUTFChars(env, model, 0);
  if (data != 0)
    chr_data = (*env)->GetStringUTFChars(env, data, 0);
  if (output != 0)
    chr_output = (*env)->GetStringUTFChars(env, output, 0);

  lp = lpx_read_model((char*)chr_model, (char*)chr_data, (char*)chr_output);

  if (model != 0)
    (*env)->ReleaseStringUTFChars(env, model, chr_model);
  if (data != 0)
    (*env)->ReleaseStringUTFChars(env, data, chr_data);
  if (output != 0)
    (*env)->ReleaseStringUTFChars(env, output, chr_output);

  if (lp == 0)
    return (0);

  mid = (*env)->GetMethodID(env, cls, "<init>", "()V");
  if (mid == 0) {
    lpx_delete_prob(lp);
    return (0);
  }
  ret = (*env)->NewObject(env, cls, mid);
  if ((*env)->ExceptionOccurred(env)) {
    lpx_delete_prob(lp);
    return (0);
  }

  put_lpx(env, ret, lp);
  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_writeMps(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_write_mps(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_writeLpt(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_write_lpt(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_printProb(JNIEnv *env,
  jobject obj, jstring fname)
{
  jint ret;

  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_print_prob(get_lpx(env, obj), (char*)chr);
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_readBas(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_read_bas(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_writeBas(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_write_bas(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_printSol(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_print_sol(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_printIps(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_print_ips(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

JNIEXPORT jboolean JNICALL
Java_org_gnu_glpk_GlpkSolver_printMip(JNIEnv *env, jobject obj, jstring fname) {
  jboolean ret;
  
  const char* chr = (*env)->GetStringUTFChars(env, fname, 0);
  ret = lpx_print_mip(get_lpx(env, obj), (char*)chr) == 0;
  (*env)->ReleaseStringUTFChars(env, fname, chr);

  return (ret);
}

/* Functions missing from documentation */

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_reallocProb(JNIEnv *env, jobject obj, jint m_max, jint n_max) {
  lpx_realloc_prob(get_lpx(env, obj), (int)m_max, (int)n_max);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_reduceForm(JNIEnv *env, jobject obj, jint len,
  jintArray ndx, jdoubleArray val, jdoubleArray work)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);
  jdouble *workEls = (*env)->GetDoubleArrayElements(env, work, 0);

  ret = lpx_reduce_form(lp, (int)len, (int*)ndxEls,
                                      (double*)valEls, (double*)workEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, work, workEls, 0);

  return ((jint)ret);
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_primOpt(JNIEnv *env, jobject obj) {
  return (jint)lpx_prim_opt(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_primArt(JNIEnv *env, jobject obj) {
  return (jint)lpx_prim_art(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_dualOpt(JNIEnv *env, jobject obj) {
  return (jint)lpx_dual_opt(get_lpx(env, obj));
}

JNIEXPORT jint JNICALL
Java_org_gnu_glpk_GlpkSolver_mixedGomory(JNIEnv *env, jobject obj, jintArray kind,
  jint len, jintArray ndx, jintArray val, jintArray work)
{
  int ret;
  LPX *lp = get_lpx(env, obj);
  jint *kindEls = (*env)->GetIntArrayElements(env, kind, 0);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);
  jdouble *workEls = (*env)->GetDoubleArrayElements(env, work, 0);

  ret = lpx_mixed_gomory(lp, (int*)kindEls, (int)len, (int*)ndxEls,
                             (double*)valEls, (double*)workEls);

  (*env)->ReleaseIntArrayElements(env, kind, kindEls, 0);
  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, work, workEls, 0);

  return ((jint)ret);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_evalRedCost(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val)
{
  double ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_eval_red_cost(lp, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return ((jdouble)ret);
}

JNIEXPORT jdouble JNICALL
Java_org_gnu_glpk_GlpkSolver_evalActivity(JNIEnv *env, jobject obj,
  jint len, jintArray ndx, jdoubleArray val)
{
  double ret;
  LPX *lp = get_lpx(env, obj);
  jint *ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  jdouble *valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  ret = lpx_eval_activity(lp, (int)len, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  return ((jdouble)ret);
}

JNIEXPORT void JNICALL
Java_org_gnu_glpk_GlpkSolver_estimObjChg(JNIEnv *env, jobject obj, jint k,
  jdouble dn_val, jdouble up_val, jintArray ndx, jintArray val, jobject chg)
{
  LPX *lp = get_lpx(env, obj);
  jclass cls = (*env)->GetObjectClass(env, chg);
  jint *ndxEls;
  jdouble *valEls;
  double dn_chg, up_chg;
  jfieldID fidDnChg, fidUpChg;
  if ((fidDnChg = (*env)->GetFieldID(env, cls, "dnChg", "D")) == 0 ||
      (fidUpChg = (*env)->GetFieldID(env, cls, "upChg", "D")) == 0)
    return;

  ndxEls = (*env)->GetIntArrayElements(env, ndx, 0);
  valEls = (*env)->GetDoubleArrayElements(env, val, 0);

  lpx_estim_obj_chg(lp, (int)k, (double)dn_val, (double)up_val,
                        &dn_chg, &up_chg, (int*)ndxEls, (double*)valEls);

  (*env)->ReleaseIntArrayElements(env, ndx, ndxEls, 0);
  (*env)->ReleaseDoubleArrayElements(env, val, valEls, 0);

  (*env)->SetDoubleField(env, chg, fidDnChg, (jdouble)dn_chg);
  (*env)->SetDoubleField(env, chg, fidUpChg, (jdouble)up_chg);
}


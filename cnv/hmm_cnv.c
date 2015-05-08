 //#include"/home/garyc/development/mpi_lasso/src/hmm_dimension.h"
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#include"hmm_cnv_dimension.h"
//#include"/home/garykche/gary/code/src/hmm_cnv_dimension.h"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void rescale_transition(
__const int states,
__constant int * cn_vec,
__constant float * forward_gap,
__constant float * backward_gap,
__global float * penalty_mat,
__global float * trans_mat,
__global float * forward_scaled_mat,
__global float * backward_scaled_mat,
__local float * rowsum,
__local float * rowvals
)
{
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int direction = get_group_id(1);
  int futurestate = get_local_id(0);

  rowsum[futurestate] = 0;
  rowvals[futurestate] = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
    float gap = direction==0?*forward_gap:*backward_gap;
  if (futurestate<states ){
   rowvals[futurestate] = rowsum[futurestate] = trans_mat[currentstate*states+futurestate] * (1-exp(-gap/HMM_TRANSITION_LARGE)) ;
   //rowvals[futurestate] = rowsum[futurestate] = trans_mat[currentstate*states+futurestate] * (penalty_mat[currentstate*states+futurestate]) * (1-exp(-gap/HMM_TRANSITION_LARGE));
   //rowvals[futurestate] = rowsum[futurestate] = trans_mat[currentstate*states+futurestate] * (1-(penalty_mat[currentstate*states+futurestate]) * exp(-gap/HMM_TRANSITION_LARGE));
  }
  barrier(CLK_LOCAL_MEM_FENCE);
   rowsum[currentstate] = 0;
  //rowvals[currentstate] = trans_mat[currentstate*states+currentstate];
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (futurestate<s) {
      rowsum[futurestate] += rowsum[futurestate+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  rowvals[currentstate] = 1-rowsum[0]>HMM_SMALL?(1-rowsum[0]):rowvals[currentstate];
  //rowvals[currentstate] *= rowsum[0]>(1-HMM_TRANSITION_DIAG)?rowsum[0]/(1-HMM_TRANSITION_DIAG):1;
  barrier(CLK_LOCAL_MEM_FENCE);
  //rowsum[0]+=rowvals[currentstate];
  if (direction==0 && futurestate<states){
    forward_scaled_mat[currentstate*states+futurestate] = rowvals[futurestate];
  }else if (direction==1 && futurestate<states){ 
    backward_scaled_mat[currentstate*states+futurestate] = rowvals[futurestate];
  }else{
     return;
  }

}

__kernel void forward_backward(
__const int states,
__const int observations,
__global int * counter,
__global float * forward_scaled_trans_mat,
__global float * backward_scaled_trans_mat,
__global float * emission_cache_mat,
__global double * forward_mat_log,
__global double * backward_mat_log,
  __local double * avg_log,
  __local double * sum_log
)
{
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int direction = get_group_id(1);
  int otherstate = get_local_id(0);
  int obs = counter[currentstate];
  avg_log[otherstate] = 0;
  sum_log[otherstate] = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
  if(direction==0 && otherstate<states){
    avg_log[otherstate] = sum_log[otherstate] = forward_mat_log[(obs-1)*states+otherstate] + log(forward_scaled_trans_mat[otherstate*states+currentstate])+log(emission_cache_mat[obs*states+currentstate]);
    barrier(CLK_LOCAL_MEM_FENCE);
  }else if(direction==1 && otherstate<states){
    avg_log[otherstate] = sum_log[otherstate] = log(backward_scaled_trans_mat[currentstate*states+otherstate])+log(emission_cache_mat[(observations-obs)*states+otherstate]) + backward_mat_log[(observations-obs)*states+otherstate];
    barrier(CLK_LOCAL_MEM_FENCE);
  }else{
    return;
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (otherstate<s) {
      avg_log[otherstate] += avg_log[otherstate+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  double avg = -avg_log[0]/(states);
  sum_log[otherstate] = exp(sum_log[otherstate]+avg);
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (otherstate<s) {
      sum_log[otherstate] += sum_log[otherstate+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (otherstate==0){
    if(direction==0){
      //forward_mat_log[obs*states+currentstate] = -log(sum_log[0]);
      forward_mat_log[obs*states+currentstate] = log(sum_log[0])-avg;
    }else if(direction==1){
      //backward_mat_log[(observations-obs-1)*states+currentstate] = -log(sum_log[0]);
      backward_mat_log[(observations-obs-1)*states+currentstate] = log(sum_log[0])-avg;
      counter[currentstate]+=1;
    }else{
      return;
    }
  }
  return;
}

__kernel void viterbi(
__const int states,
__const int observations,
__global int * counter,
__global float * forward_scaled_trans_mat,
__global float * emission_cache_mat,
__global double * greedy_mat_log,
__global int * bestpath_mat,
  __local int * best_index,
  __local double * best_prob
)
{
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int otherstate = get_local_id(0);
  int obs = counter[currentstate];

  best_index[otherstate] = 0;
  best_prob[otherstate] = -1e11;
  barrier(CLK_LOCAL_MEM_FENCE);
  if(otherstate<states){
    best_prob[otherstate] = greedy_mat_log[(obs-1)*states+otherstate] + log(forward_scaled_trans_mat[otherstate*states+currentstate])+log(emission_cache_mat[obs*states+currentstate]);
    best_index[otherstate] = otherstate;
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (otherstate<s && best_prob[otherstate]<best_prob[otherstate+s]){
      best_prob[otherstate] = best_prob[otherstate+s];
      best_index[otherstate] = best_index[otherstate+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (otherstate==0){
    greedy_mat_log[obs*states+currentstate] = best_prob[0];
    bestpath_mat[obs*states+currentstate] = best_index[0];
    counter[currentstate]+=1;
  }
  return;
}

__kernel void likelihood(
__const int states,
__global double * forward_mat_log,
__global double * backward_mat_log,
__global float * likelihood_log,
__local double * temp1,
__local double * temp2,
__local double * temp3
)
{
  int state = get_local_id(0);
  temp1[state] = 0;
  temp2[state] = 0;
  temp2[state] = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
  if (state<states){
    temp1[state] = temp2[state] = forward_mat_log[state] + backward_mat_log[state];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  // compute average log value for rescaling
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (state<s) {
      temp1[state] += temp1[state+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  double avg;
  //avg = -temp1[0]/9;
  avg = -temp1[0]/(states);
  // *likelihood_log = avg;
  //return;
  if(state<states) temp3[state] = exp(temp2[state]+avg);
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int s=HMM_STATE_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (state<s) {
      temp3[state] += temp3[state+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  *likelihood_log = log(temp3[0])-avg;
  return;
}

__kernel void posterior(
__const int states,
__global double * forward_mat_log,
__global double * backward_mat_log,
__constant float * likelihood_log,
__global double * prob_mat_log,
__local double * temp1,
__local double * temp2
)
{
  int obs = get_group_id(0);
  int state = get_local_id(0);
  if (state<states){
    temp1[state] = forward_mat_log[obs*states+state] ;
    temp2[state] = backward_mat_log[obs*states+state];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (state<states) prob_mat_log[obs*states+state] = temp1[state]+temp2[state]-*likelihood_log;
}

__kernel void baumwelch_normalize(
__const int states,
__const int observations,
__global double * prob_mat_log,
__global double * baumwelch_norm,
__local double * temp1
)
{
  int state = get_group_id(0);
  if (state>=states) return;
  int tid = get_local_id(0);
  temp1[tid] = 0;
  barrier(CLK_LOCAL_MEM_FENCE);
  //for(int chunk=0;chunk<2;++chunk){    
  for(int chunk=0;chunk<(observations-1)/HMM_OBS_BLOCK_WIDTH+1;++chunk){    
    int index = chunk*HMM_OBS_BLOCK_WIDTH+tid;
    if (index<(observations-1)){
      temp1[tid] += exp(prob_mat_log[index*states+state]);
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  } 
  if(tid==-1){
    float sum = 0;
    for(int i=0;i<HMM_OBS_BLOCK_WIDTH;++i) sum+=temp1[i];
    baumwelch_norm[state] = sum;
    return;
  }
  for(int s=HMM_OBS_BLOCK_WIDTH/2; s>0; s>>=1) {
    if (tid<s) {
      temp1[tid] += temp1[tid+s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if(tid==0){
    baumwelch_norm[state] = temp1[0];
  }
  return;
}

__kernel void baumwelch_transition(
__const int states,
__global int * counter,
__global double * forward_mat_log,
__global double * backward_mat_log,
__global float * emission_cache_mat,
__constant float * likelihood_log,
__global float * running_total_trans_mat,
__global float * scaled_trans_mat
)
{
  int obs = *counter;
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int futurestate = get_local_id(0);
  if (futurestate<states){
    running_total_trans_mat[currentstate*states+futurestate] += exp(forward_mat_log[obs*states+currentstate]+log(scaled_trans_mat[currentstate*states+futurestate])+log(emission_cache_mat[(obs+1)* states+futurestate]) + backward_mat_log[(obs+1)*states+futurestate] - *likelihood_log);
  }
  *counter = obs+1;
}

__kernel void zero_transition(
__const int states,
__global float * new_trans_mat
)
{
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int futurestate = get_local_id(0);
  if (futurestate<states) new_trans_mat[currentstate*states+futurestate] = 0;
}

__kernel void normalize_transition(
__const int states,
__global float * trans_mat,
__global float * new_trans_mat,
__constant double * baumwelch_norm
)
{
  int currentstate = get_group_id(0);
  if (currentstate>=states) return;
  int futurestate = get_local_id(0);
  float t = futurestate<states?new_trans_mat[currentstate*states+futurestate]/baumwelch_norm[currentstate]:1;
  barrier(CLK_LOCAL_MEM_FENCE);
  t = HMM_SMALL + (1-HMM_SMALL) * t;
  if (futurestate<states) trans_mat[currentstate*states+futurestate] = t;
}


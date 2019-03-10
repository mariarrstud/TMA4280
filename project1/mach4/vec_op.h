#ifndef VEC_OP_H_
#define VEC_OP_H_

void print_vec(double *vec, int n);
double sum_vec_using_openmp(double *vec, int n);
void slice_vec(double *vec, double *part_vec, int start, int end);

#endif /* VEC_OP_H_ */

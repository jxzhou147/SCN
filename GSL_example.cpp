#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/*����ODE�Ҳ෽����ʽ*/
int  func (double t, const double y[], double f[], void *params)
{
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

/*����ODE�ſɱȾ�����ʽ*/

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
   double mu = *(double *)params;
   /*������ε���˼��*/
   gsl_matrix_view dfdy_mat= gsl_matrix_view_array (dfdy, 2, 2);  //������ͼ 
   gsl_matrix * m = &dfdy_mat.matrix;   //ָ�룬ȡ����Ԫ��ַ�õ� 
   gsl_matrix_set (m, 0, 0, 0.0);    //���ʾ����еĵ�i+1�У�j+1�е�Ԫ��:gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);
   gsl_matrix_set (m, 0, 1, 1.0);
   gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
   gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
   dfdt[0] = 0.0;
   dfdt[1] = 0.0;
   return GSL_SUCCESS;
}


//�� evolution function�任�������� Van der Pol oscillator
  int main (void)
{
	
   const gsl_odeiv_step_type * T= gsl_odeiv_step_rk8pd;    //����ά�������8�׷��� 
   gsl_odeiv_step * s= gsl_odeiv_step_alloc (T, 2);    //����һ��2άϵͳ������ΪT��stepping���� ��ÿ���ظ�ʹ�ý�������sʱ���á� 
   gsl_odeiv_control * c= gsl_odeiv_control_y_new (1e-6, 0.0); //�ú��������µĿ���Ŀ�꣬��Ŀ��άÿһ����ľֲ����С�ھ������������� 
   gsl_odeiv_evolve * e= gsl_odeiv_evolve_alloc (2);   //�ú�������һ��2άϵͳ��ָ���·���Ľ��������Ķ��� 
   
   
   double mu = 10;
   gsl_odeiv_system sys = {func, jac, 2, &mu};  //�������� gsl_odeiv_systemָ��壬�����������ͣ�func��jac,������ά����ָ��ϵͳ���������ָ�룩 
   
   /*����ʱ�䷶Χ��������ֵ�趨*/ 
   double t = 0.0, t1 = 100.0;
   double h = 1e-6;       //���� 
   double y[2] = { 1.0, 0.0 };  //��ֵ 
   
   FILE *ff=fopen("trace_with_jac.txt","a+");
   
while (t < t1)
{
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);        /*����ָ��*/   //�ú����ý������������ӳ�ֵ��ϵͳ�ģ�e,dydt���ƽ� 
     if (status != GSL_SUCCESS)
     break;     /*����ʧ�����Ƴ�*/ 
     printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);   /*����ɹ���������*/ 
     fprintf(ff,"%.5e %.5e %.5e\n", t, y[0], y[1]);
}
   fclose(ff);

  /*�ͷ��ڴ�*/ 
   gsl_odeiv_evolve_free (e);    //�ͷŹ����ݻ�����e�����м��� 
   gsl_odeiv_control_free (c);   //�ͷŹ��ڿ��ƺ��� c�����м��� 
   gsl_odeiv_step_free (s);      //�ͷŹ��ڽ�������s�����м��� 
   return 0;
}

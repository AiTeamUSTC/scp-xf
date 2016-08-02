#include<stdio.h> 
#include<stdlib.h>
#include<fstream>
#include<memory.h>
#include<vector>
#include<algorithm>
#include<iostream>
#include<limits.h>
#include<float.h>
//#include<unistd.h>
#include<time.h> 
using namespace std;
#define pop_size 4
#define MAXC 50000
#define MAXR 500000
#define max_step 1000

/********�������Ƭ��һÿ������λ�Ľṹ�壬ÿ������λ�����ܶ���Ϣ�����Ը��ǵ��У����еĴ��ۣ�ʱ����ȵ�********/
typedef struct gene{
	int config;  //�Ƿ���Լ����ѡ�� 
	int flag;    //�Ƿ�����Ƴ���ѡ�� 
	int score; //�еķ�������ʼֵΪ�����ܸ��ǵ��еĸ������ڵ�λ���ۼ��ϸ��������У�score��ֵԽ����У�Խ���ȱ�ѡ���ѡ��
	int time_stamp; //ʱ��� 
	int is_in_c;    // �����ѡ�������Ϊ1������Ϊ0 
	int cost;  //ÿһ�еĴ��ۣ�ʵ������ǰ2000�����ݣ���ʾ�ľ���ÿһ�еĴ��ۣ���λ���ۼ��ϸ��������У�ֵͳһΪ1
	vector<int> cover_rs;// ��Ÿ��п��Ը��ǵ���
}gene;


/********�������Ƭ�ζ�ÿ������λ�Ľṹ�壬ÿ������λ����һ�У��������е�Ȩֵ�����Լ�ÿһ�п��Ա���Щ�и���********/
typedef struct meme{
	int weight;  // ÿһ�е�Ȩ�أ����������в��ϱ仯��Խ�Ѹ��ǵ����У�ȨֵԽ�󣬴Ӷ������䱻���ǵĿ�����
	vector<int> covered_cs; // ��ſ��Ը��Ǹ��е���
}meme;


/********������Ⱥ��ÿ������Ľṹ��********/
typedef struct individual{
	gene sol_gene[MAXC];  //����Ƭ��һ
	meme sol_meme[MAXR]; //����Ƭ�ζ�
	double prior[MAXC];       //�����������Ƭ��һ����Ӧ����λ��ѡ���ѡ�������Ȩ,��λ���ۼ��ϸ��������У���score����һ��
	vector<int> uncover_rs;   //���ÿ����û�и��ǵ�����
	int fitness;                       // �򵥱�ʾ��ѡ��Ĵ�С
}individual;


/********����ȫ�ֱ���************/
int row_num;//����
int col_num;//����
double Cr = 0.8; // ������
individual ind[pop_size];//ԭ��Ⱥ
individual Mind[pop_size];//ͻ����Ⱥ
/*********�����Ķ���***********************/
double random(double low, double high);
void build_instance(char* file);
void initialize(individual * ind);
void change_in(int j, individual *Ind);
void chang_out(int j, individual*Ind);
int evaluate(individual *ind);
int find_maxc_out(individual *ind);
int find_maxc_in(individual *ind, int r);
void crossover();
void mutation();
void SLS(individual *ind);
bool cmp(individual a, individual b);


/*********������*******************/
int main(int argc, char *argv[])
{
	 //srand((unsigned)time(NULL));//��������� 
	char filename[] = "scp41.txt";
	//��Ⱥ��ģ 
	build_instance(filename);
	srand(40);
	for (int i = 0; i < pop_size; i++)
	{
		//��Ⱥ��ʼ�� 
		  initialize(&ind[i]);
		//������ʼ��Ⱥ 
		ind[i].fitness = evaluate(&ind[i]);
		cout << ind[i].fitness << endl;
	}
	for (int gen = 0; gen < max_step; gen++)
	{
		//cout << "MDZZ" << endl;
		
		mutation();
		cout << "gen: " << gen << endl;
		for (int i = 0; i < pop_size; i++)
			SLS(&ind[i]);
		crossover();
		for (int i = 0; i < pop_size; i++)
		{
			ind[i].fitness = evaluate(&ind[i]);
			cout << ind[i].fitness << endl;
		}
	}
	//sort(ind, ind+pop_size, cmp);
	for (int i = 0; i < pop_size; i++)
		cout << ind[i].fitness << endl;
	return 0;
}
/*********���������*******************/


/****��������� ***/
double random(double low, double high)
{
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}

/******ʵ����ģ********/
void build_instance(char* file)
{
	for (int it = 0; it < pop_size; it++)
	{
		int i, j, k, t, p, cm, temp;
		ifstream f1;
		f1.open(file);
		f1 >> row_num >> col_num;		//�����ܵ������Լ����� 
		
		/**************��ʼ���׶�*************************/ 
		ind[it].uncover_rs.clear();		//��ʼ���uncover_rs 
		for (i = 1; i <= row_num; i++)
		{
			ind[it].sol_meme[i].covered_cs.clear();	         //��ʼ��������е�covered_cs
		}
		for (j = 1; j <= col_num; j++)
		{
			f1 >> k;		//����ÿһ�еĴ��� 
			/*�����Ҫ�����λ���۵ļ��ϸ������⣬ֱ�Ӱ�K��ֵΪ1����*/
			//ind[it].sol_gene[j].cost = k;
			ind[it].sol_gene[j].cost = 1;
			ind[it].sol_gene[j].config = 1;    //��ʾ�����Լ����ѡ�� 
			ind[it].sol_gene[j].flag = 1;     // ��ʾ��ʼ״̬�����ж����ԴӺ�ѡ���Ƴ� 
			ind[it].sol_gene[j].score = 0;    // ��ʼ��scoreΪ0 
			ind[it].sol_gene[j].is_in_c = 0;   //��ʼ�����ж����ں�ѡ���� 
			ind[it].sol_gene[j].time_stamp = 0; // ʱ���
			ind[it].sol_gene[j].cover_rs.clear(); 			//��ʼ���cover_rs 
		}
		ind[it].fitness = 0; 		// ��ʼ���������Ӧֵ���� 
		
	   /****************��ʼ������Ӧ�����й�ϵ********************/ 
		for (i = 1; i <= row_num; i++)
		{
			ind[it].uncover_rs.push_back(i);			//��ʼ�����е���Ϊδ������
			ind[it].sol_meme[i].weight = 1;			//��ʼ��ÿһ�е�Ȩ�� 
			f1 >> t;			//�����i�п��Ա������и��ǣ�������t��  
			if (t == 1)// �������ĳһ��ֻ�ܱ�һ�и��ǣ��򽫸��й̶��ں�ѡ���� 
			{
				ind[it].sol_gene[cm].flag = 0;   // ��ʾ���в����Ա��Ƴ���ѡ��
				ind[it].sol_gene[cm].is_in_c = 1; // ��ʾ�����ں�ѡ���� 
				ind[it].uncover_rs.pop_back(); //��Ӧ�����ǵ���һ���Ƴ�uncover_rs
			}
			for (p = 0; p < t; p++)
			{
				f1 >> cm;			//������Ը��ǵ�i�е���
				ind[it].sol_meme[i].covered_cs.push_back(cm);				//��covered_cs�д���cm�� 
				ind[it].sol_gene[cm].cover_rs.push_back(i);				//��cm�е�cover_rs�д������ 
			}
		}
		//���������е�Ȩֵ�Լ�����Ȩ 
		for (j = 1; j <= col_num; j++)
		{
			for (i = 0; i < ind[it].sol_gene[j].cover_rs.size(); i++)
			{
				temp = ind[it].sol_gene[j].cover_rs[i];
				ind[it].sol_gene[j].score += ind[it].sol_meme[temp].weight;
			}
			ind[it].prior[j] = ind[it].sol_gene[j].score * 1.0 / ind[it].sol_gene[j].cost;
			//cout<<ind[it].prior[j]<<" ";
		}
		f1.close();
	}
		/*for(int i = 1;i <= row_num; i++) 
		{
			cout<<ind[0].sol_meme[i].covered_cs.size()<<endl;
			for(int j = 0; j<ind[0].sol_meme[i].covered_cs.size(); j++)
			cout<<ind[0].sol_meme[i].covered_cs[j]<<" ";
			cout<<endl;
		}	*/
}

/*****�γɳ�ʼ�� ******/
void initialize(individual * ind)
{
		   int best_array[MAXC]; 
		   double max;
		   int cnt; 
		   int i, j;
		   while (!(*ind).uncover_rs.empty())
		   {
				   //srand((unsigned)time(NULL));//���������

			       max = -DBL_MAX;
			       cnt = 0;
			       for (j = 1; j <= col_num; j++)
				   {
					   if ((*ind).sol_gene[j].is_in_c)
							  continue;
					   if (max<(*ind).prior[j])
				        {
						   max = (*ind).prior[j];
					           best_array[0] = j;
					           cnt = 1;
				        }
					   else if (max == (*ind).prior[j])
					           best_array[cnt++] =j;
			        }
			        if (cnt>0)
				    {
				        vector<int> bestA;
				        int minCost = INT_MAX;
				        for (i = 0; i < cnt; i++)
						{
					        if (best_array[i] < minCost)
							{
						        minCost = best_array[i];
						        bestA.clear();
						        bestA.push_back(i);
					        }  
					          else if (minCost == best_array[i]) 
						        bestA.push_back(i);
				        }
			            i= rand() % bestA.size();
			            //cout<<i<<" ";
			            cout<<best_array[bestA[i]]<<" ";
				        change_in(best_array[bestA[i]], ind);
			        }
		    }
		    cout<<endl;
}
/*change_in��������Ҫ�Ի���Ƭ��һ�л���λ��is_in_c�������ı䣬is_in_c��0��1��ʾ����Ӧ�м����ѡ�⣬������Ҳ�������ھӵ�����Ȩ����*/ 
void change_in(int j, individual  *Ind)
{
	   int temp_r, temp_c, btemp_c, temp, count;
 	   (*Ind).sol_gene[j].is_in_c = 1;	   //��Ӧ����λ���Ա���Ϊ1,��ʾ��j�м����ѡ��
	   (*Ind).sol_gene[j].score = -(*Ind).sol_gene[j].score; 	   // score��Ϊ�෴��  
	   (*Ind).prior[j] = -(*Ind).prior[j]; 	   //��Ӧ�ĸ�������Ȩ 
	   	/*�ҳ���j�е��ھӣ��������ھӵ�����Ȩ�Լ�״̬*/ 
       for (int i = 0; i < (*Ind).sol_gene[j].cover_rs.size(); i++)
	   {
		   temp_r = (*Ind).sol_gene[j].cover_rs[i];     //��j�п��Ը��ǵ��� 
		   count = 0;  			//���ڼ���Ƿ��ǵ�һ�θ���temp_r 
		   /*�ҵ�j�й���temp_r���ھ� �� ����Ȩֵ��״̬ */
		   for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
		   {
			   temp_c = (*Ind).sol_meme[temp_r].covered_cs[k]; 				 //���Ը���temp_r�е��� 
			   if (temp_c == j)   continue;				 //�����j��ֱ������ 
			   /*�ж�temp_c���Ƿ��ں�ѡ����*/
			   if ((*Ind).sol_gene[temp_c].is_in_c == 1)
			   {
				   btemp_c = temp_c; //��¼���ں�ѡ���е��ھ� 
				   count++;
			   }

			   (*Ind).sol_gene[temp_c].config = 1;  		//���������ھ�״̬
		   }
		   /*countΪ0��ʾ��temp_r��һ�θ��ǣ���Ӧ�ĵ�j�й���temp_r�������ھӵ�score���ȥtemp_r�е�Ȩֵ������������Ȩ*/ 
		   if (count == 0)
		   {
			   for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
			   {
				   temp = (*Ind).sol_meme[temp_r].covered_cs[k];
				   if (temp == j)
					   continue;
				   (*Ind).sol_gene[temp].score -= (*Ind).sol_meme[temp_r].weight;
				   (*Ind).prior[temp] = (*Ind).sol_gene[temp].score * 1.0 / (*Ind).sol_gene[temp].cost;
			   }

		   }
		   /*�ں�ѡ���е��ھ�score��Ӧ���Ǹ�ֵ������ÿ��һ�����ھӽ����ѡ�⣬��Ӧ��temp_r��Ȩֵ��Ӱ�����Ҫ��Ĩȥ*/
		   /*����ΪʲôֻĨȥ��һ�������Ľ����ǵ�temp_r����θ���ʱ��temp_r��Ȩֵ����score��Ӱ����Բ���*/ 
		   else if (count == 1)
		   {
			   (*Ind).sol_gene[btemp_c].score += (*Ind).sol_meme[temp_r].weight;
			   (*Ind).prior[btemp_c] = (*Ind).sol_gene[btemp_c].score * 1.0 / (*Ind).sol_gene[btemp_c].cost;
		   }
		   //��temp_r��uncover_rs��ɾ��
		   vector<int>::iterator it;
		   it = (*Ind).uncover_rs.begin();
		   while (it != (*Ind).uncover_rs.end())
		   {
			   if (*it == temp_r)
			   {
				   it = (*Ind).uncover_rs.erase(it);
				   break;
			   }
			   else
				   it++;
		   }
	   }
}
/*change_out��������Ҫ�Ի���Ƭ��һ�л���λ��is_in_c�������ı䣬is_in_c��1��0��ʾ����Ӧ���Ƴ���ѡ�⣬������Ҳ�������ھӵ�����Ȩ����*/ 
void change_out(int j, individual * ind)
{
	int temp_r;
	int temp_c;
	bool first;
	//��Ӧ����λ����Ϊ0
	(*ind).sol_gene[j].is_in_c = 0;
	//����score 
	(*ind).sol_gene[j].score = -(*ind).sol_gene[j].score;
	//��������Ȩ 
	(*ind).prior[j] = -(*ind).prior[j];
	(*ind).sol_gene[j].config = 0;
	//�����ھӵ�״̬, ��ʾ���Ա������ѡ�� 
	for (int i = 0; i < (*ind).sol_gene[j].cover_rs.size(); i++)
	{
		//��j�п��Ը��ǵ��� 
		temp_r = (*ind).sol_gene[j].cover_rs[i];
		first = true;
		//�ҵ�j�й���temp_r���ھ� 
		for (int k = 0; k < (*ind).sol_meme[temp_r].covered_cs.size(); k++)
		{
			temp_c = (*ind).sol_meme[temp_r].covered_cs[k];
			if (temp_c == j) continue;
			else
			{
				if ((*ind).sol_gene[temp_c].is_in_c == 1)
				{
					first = false;
				}
			}
			(*ind).sol_gene[temp_c].config = 1;
		}
		//����ѡ������Ȼ���ڿ��Ը���temp_r����ʱ 
		if (!first)
		{
			//�ҳ����Ը���temp_r�������� 
			for (int k = 0; k < (*ind).sol_meme[temp_r].covered_cs.size(); k++)
			{
				temp_c = (*ind).sol_meme[temp_r].covered_cs[k];
				if (temp_c == j) continue;
				//���º�ѡ������ھӵ�Ȩֵ 
				else if (temp_c != j && (*ind).sol_gene[temp_c].is_in_c == 0)
				{
					(*ind).sol_gene[temp_c].score -= (*ind).sol_meme[temp_r].weight;
					(*ind).prior[temp_c] = (*ind).sol_gene[temp_c].score * 1.0 / (*ind).sol_gene[temp_c].cost;
				}
				//���º�ѡ���ڵ��ھӵ�Ȩֵ 
				else if (temp_c != j && (*ind).sol_gene[temp_c].is_in_c == 1)
				{
					(*ind).sol_gene[temp_c].score += (*ind).sol_meme[temp_r].weight;
					(*ind).prior[temp_c] = (*ind).sol_gene[temp_c].score * 1.0 / (*ind).sol_gene[temp_c].cost;
				}
			}
		}
		else if (first)
		{
			(*ind).uncover_rs.push_back(temp_r);
			//�����ھӵ�Ȩֵ 
			for (int k = 0; k < (*ind).sol_meme[temp_r].covered_cs.size(); k++)
			{
				temp_c = (*ind).sol_meme[temp_r].covered_cs[k];
				if (temp_c == j) continue;
				else
				{
					(*ind).sol_gene[temp_c].score += (*ind).sol_meme[temp_r].weight;
					(*ind).prior[temp_c] = (*ind).sol_gene[temp_c].score * 1.0 / (*ind).sol_gene[temp_c].cost;
				}
			}
		}
	}
}
/*ͻ���������Ҫ���γ�ͻ����Ⱥ��ͻ�䷢���ڻ���Ƭ�ζ��������������ͻ����Ⱥ��ԭ��Ⱥ��������Ҫ�ο���һ��DE�㷨��ͻ���������*/ 
void mutation()
{

	//ͨ��ͻ��Ĳ��������µ���Ⱥ 
	for (int i = 0; i < pop_size; i++)
	{
		for (int j = 1; j <= col_num; j++)
		{

			Mind[i].sol_gene[j].config = ind[i].sol_gene[j].config;
			Mind[i].sol_gene[j].cost = ind[i].sol_gene[j].cost;
			Mind[i].sol_gene[j].time_stamp = ind[i].sol_gene[j].time_stamp;
			Mind[i].sol_gene[j].flag = ind[i].sol_gene[j].flag;
			Mind[i].sol_gene[j].is_in_c = ind[i].sol_gene[j].is_in_c;
			Mind[i].sol_gene[j].score = ind[i].sol_gene[j].score;
			Mind[i].prior[j] = ind[i].prior[j];
			Mind[i].sol_gene[j].cover_rs.clear();
			for (int s = 0; s < ind[i].sol_gene[j].cover_rs.size(); s++)
			{
				Mind[i].sol_gene[j].cover_rs.push_back(ind[i].sol_gene[j].cover_rs[s]);
			}
		}
		for (int ri = 1; ri <= row_num; ri++)
		{

			Mind[i].sol_meme[ri].weight = ind[i].sol_meme[ri].weight;
			for (int s = 0; s < ind[i].sol_meme[ri].covered_cs.size(); s++)
			{
				Mind[i].sol_meme[ri].covered_cs.push_back(ind[i].sol_meme[ri].covered_cs[s]);
			}
		}
		for (int s = 0; s <ind[i].uncover_rs.size(); s++)
		{
			Mind[i].uncover_rs.push_back(ind[i].uncover_rs[s]);
		}
		Mind[i].fitness = ind[i].fitness;
	}
	//ͨ��ͻ����������µ���Ⱥ 
	int r1, r2, r3;
	for (int k = 0; k < pop_size; k++)
	{
		int i = rand() % pop_size;
		do
		{
			r1 = rand() % pop_size;
		} while (r1 == i);
		do
		{
			r2 = rand() % pop_size;
		} while (r2 == i || r2 == r1);
		do
		{
			r3 = rand() % pop_size;
		} while (r3 == i || r3 == r1 || r3 == r2);
		for (int j = 1; j <= col_num; j++)
		{
			int temp;
			temp = 0;
			temp = Mind[r1].sol_gene[j].is_in_c + Mind[r2].sol_gene[j].is_in_c + Mind[r3].sol_gene[j].is_in_c;
			if (temp <= 1 && Mind[i].sol_gene[j].flag == 1)
			{
				change_out(j, &Mind[i]);
			}
			else if (temp >= 2 && Mind[i].sol_gene[j].config == 1)
			{
				change_in(j, &Mind[i]);
			}
		}
	}
}
/*��������������ind�еĽ��ĳЩ�У�ִ�н������ӣ�������meme����Ƭ�ν��н��棬ֱ�ӽ����������indȺ����*/
/*����ͻ�����������������һ����Ⱥ����Ϣ���������ã�ѡ���Լ��Խ�ĸĽ���������SLS����*/ 
void crossover()
{
	for (int i = 0; i < pop_size; i++)
	{
		int j_rand = rand() % row_num;
		for (int j = 1; j <= row_num; j++)
		{
			double rand1 = random(0, 1);
			if (rand1 < Cr || j == j_rand)
			{
				ind[i].sol_meme[j].weight = Mind[i].sol_meme[j].weight;
			}
		}
	}
}
 /*��һ�������Ľ⿪ʼ��ִ������ֲ��ѹ���
  * individual* ind, һ��ָ���Ľ����
  *����ֲ������Ĳ������ð�����meme����Ƭ���У�����sls()�ĳ�������ʵ����Խ׶�Ϊ����ʱ�䣬ֻ��ȡ1000��
  *��һ���Ľ��������ʱ�����õĸ���һЩ ��RWLS��30000000 
  */
void SLS(individual *ind)
{

	int  step = 1;
	int  temp_cost = evaluate(ind);
	while (step <= 1000)
	{
	//	cout << "step: " << step << endl;
		if ((*ind).uncover_rs.empty())
		{
			int maxc;
			maxc = find_maxc_out(ind);
			//cout << "MDZZ" << endl;
			change_out(maxc, ind);
			(*ind).sol_gene[maxc].time_stamp = step;
			continue;//����ط�����while��Ϊ�˱�֤��ʼ�׶�change_out���е�time_stamp��ͬ 
		}
		int  Maxc;
		Maxc = find_maxc_out(ind);
		change_out(Maxc, ind);
		(*ind).sol_gene[Maxc].time_stamp = step;
		while (!(*ind).uncover_rs.empty())
		{
	    	int i = rand() % (*ind).uncover_rs.size();
		    int r = (*ind).uncover_rs[i];
		    int maxc = find_maxc_in(ind, r);
		    if ((*ind).sol_gene[maxc].cost + evaluate(ind) > temp_cost)
	    	break;
		    change_in(maxc, ind);
		    (*ind).sol_gene[maxc].time_stamp = step;
		    //δ�����ǵ��е�Ȩֵ���ӣ���Ӧ�Ŀ��Ը��Ǹ��е��е�score������Ȩ����
		   for (int i = 0; i < (*ind).uncover_rs.size(); i++)
		   {
		     int ur = (*ind).uncover_rs[i];
		     (*ind).sol_meme[ur].weight += 1;
		     for (int j = 0; j < (*ind).sol_meme[ur].covered_cs.size(); j++)
		     {
	               	//���Ը��Ǹ��е���
	             	int temp_j;
		            temp_j = (*ind).sol_meme[ur].covered_cs[j];
		            (*ind).sol_gene[temp_j].score += 1;
		            (*ind).prior[temp_j] = (*ind).sol_gene[temp_j].score * 1.0 / (*ind).sol_gene[temp_j].cost;
		     }
	    	}
		}
		step++;
	}
}
/**
  *brief �����/�������Ӧ��ֵ
  *param[in] individual* ind, һ��ָ���Ľ����
  *return int fitness, ��ѡ�е��е���Ŀ
  */
int evaluate(individual *ind)
{
	int fitness;
	fitness = 0;
	for (int j = 1; j <= col_num; j++)
	{
		if ((*ind).sol_gene[j].is_in_c == 1)
		{
			fitness += (*ind).sol_gene[j].cost;
		}
	}
	return fitness;
}
/**
  *brief �Ӻ�ѡ����ѡ��score�������Ƴ� 
  *param[in] individual* ind, һ��ָ���Ľ����
  *return maxc����ѡ���б�ѡ���Ƴ�����
  */
int find_maxc_out(individual * ind)
{
	int maxc;
	maxc = 1;
	double  max = -DBL_MAX;
	for (int j = 1; j <= col_num; j++)
	{
		if ((*ind).sol_gene[j].is_in_c == 0)
			continue;
		if ((*ind).sol_gene[j].flag == 0)
			continue;
		if (max<(*ind).prior[j])
		{
			max = (*ind).prior[j];
			maxc = j;
		}
		else if (max == (*ind).prior[j])
		{
			if ((*ind).sol_gene[j].time_stamp < (*ind).sol_gene[maxc].time_stamp)
			{
				maxc = j;
			}
		}
	}
	return maxc;
}
/**
  *brief �ӷǺ�ѡ�������ѡ��score�����������ѡ�� 
  *param[in] individual* ind, һ��ָ���Ľ����
  *return maxc����ѡ�������ѡ�����
  */
int find_maxc_in(individual* ind, int r)
{

	double  max = -DBL_MAX;
	int maxc;
	maxc = 1;
	for (int j = 0; j < (*ind).sol_meme[r].covered_cs.size(); j++)
	{
		int temp_j;
		temp_j = (*ind).sol_meme[r].covered_cs[j];
		if ((*ind).sol_gene[temp_j].is_in_c == 1)
			continue;
		if ((*ind).sol_gene[temp_j].config == 0)
			continue;
		if (max<(*ind).prior[temp_j])
		{
			max = (*ind).prior[temp_j];
			maxc = temp_j;
		}
		else if (max == (*ind).prior[temp_j])
		{
			if ((*ind).sol_gene[temp_j].time_stamp < (*ind).sol_gene[maxc].time_stamp)
			{
				maxc = temp_j;
			}
		}
	}
	return maxc;
}
/*д�������Ŀ��������Ⱥ�Ÿ���ѡ�����ŵĸ�����뵽ͻ�������У�Ŀǰ��û�õ�*/ 
bool cmp(individual  *a, individual *b)
{
	return  (*a).fitness <  (*b).fitness;
}



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
#define pop_size 20
#define MAXC 50000
#define MAXR 500000
#define max_step 100


/**
  *brief �������Ƭ��һÿ������λ�Ľṹ�壬ÿ������λ����һ�����С�/���Ӽ��ϡ������ܸ������ɡ��С�/��Ԫ�ء�
*/
typedef struct gene{
	int config;             ///�Ƿ���Լ����ѡ��
	int flag;               ///�Ƿ�����Ƴ���ѡ��
	int score;              ///�еķ�����ֵԽ���ʾ......��ֵԽС��ʾ......
	int time_stamp;         ///ʱ���
	int is_in_c;            ///�����ѡ�������Ϊ1������Ϊ0
	int cost;               ///ʲô��˼���������
	vector<int> cover_rs;   /// rows it can cover
}gene;



/**
  *brief �������Ƭ�ζ�ÿ������λ�Ľṹ�壬һ������λ����һ����Ҫ���ǡ����У�����
*/
typedef struct meme{
	int weight;             /// ��ʾһ��Ԫ�ػ��б����ǵġ��Ѷȡ���ֵԽ�󣬱�ʾ����Խ����
	vector<int> covered_cs; /// cols have can cover this row
}meme;


/**
  *brief ������Ⱥ��ÿ������Ľṹ��
*/
typedef struct individual{
	gene sol_gene[MAXC];    ///����Ƭ��1
	meme sol_meme[MAXR];    ///����Ƭ��2
	double prior[MAXC];     ///ʲô��
	vector<int> uncover_rs; ///δ�����ǵ���
	int fitness;            ///��Ӧ�ȵ�ֵ����ô���㣬��ʾʲô����ѡ�е��е���Ŀ��

}individual;



int row_num;                        ///instance ����������
int col_num;                        ///instance ����������
double Cr = 0.8;                    ///�㷨������ ������ʣ�
individual ind[pop_size];           ///�㷨������Ⱥ��
individual Mind[pop_size];          ///�㷨����������Ⱥ�壿
vector<int> best_sol;               ///�㷨ִ�й������ҵ������Ž�


///����������
double random(double low, double high);
void build_instance(char* file);
void initialize();
void change_in(int j, individual *Ind);
void chang_out(int j, individual*Ind);
int evaluate(individual *ind);
int find_maxc_out(individual *ind);
int find_maxc_in(individual *ind, int r);
void crossover();
void mutation();
void SLS(individual *ind);
bool cmp(individual a, individual b);


/**
  * @brief ������
  * @author xufang
  * @param [in] NULL
  * @param [out] NULL
  * @return NULL
  * @note ....todo
  */
int main(int argc, char *argv[])
{
	srand((unsigned)time(NULL));//���������
	char filename[] = "scp41.txt";
	//��Ⱥ��ģ
	build_instance(filename);
	//��Ⱥ��ʼ��
	initialize();
	for (int i = 0; i < pop_size; i++)
	{
		//������ʼ��Ⱥ
		ind[i].fitness = evaluate(&ind[i]);
		//cout << ind[i].fitness << endl;
	}
	for (int gen = 0; gen < max_step; gen++)
	{
		mutation();
		for (int i = 0; i < pop_size; i++)
			SLS(&ind[i]);
		crossover();
		for (int i = 0; i < pop_size; i++)
		{
			ind[i].fitness = evaluate(&ind[i]);
			cout << ind[i].fitness << endl;
		}
	}
	/*sort(ind, ind+pop_size, cmp);
	for (int i = 0; i < pop_size; i++)
		cout << ind[i].fitness << endl;*/
	return 0;
}
///���������


/**
  *brief ���������
  */
double random(double low, double high)
{
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}


/**
  *brief ʵ����ģ
  */
void build_instance(char* file)
{
	for (int it = 0; it < pop_size; it++)
	{
		int i, j, k, t, p, cm, temp;
		ifstream f1;
		f1.open(file);
		//�����ܵ������Լ�����
		f1 >> row_num >> col_num;
		//��ʼ���uncover_rs
		ind[it].uncover_rs.clear();
		//����ÿһ�еĴ���
		for (j = 1; j <= col_num; j++)
		{
			f1 >> k;
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
		// ��ʼ���������Ӧֵ����
		ind[it].fitness = 0;
		for (i = 1; i <= row_num; i++)
		{
			//��ʼ���covered_cs
			ind[it].sol_meme[i].covered_cs.clear();
			//��ʼ�����е���Ϊδ������
			ind[it].uncover_rs.push_back(i);
		}
		//��������еĸ��ǹ�ϵ
		for (i = 1; i <= row_num; i++)
		{
			//��ʼ��ÿһ�е�Ȩ��
			ind[it].sol_meme[i].weight = 1;
			//�����i�п��Ա������и���
			f1 >> t;
			//������Ը��ǵ�i�е���
			for (p = 0; p < t; p++)
			{
				f1 >> cm;
				//��covered_cs�д���cm��
				ind[it].sol_meme[i].covered_cs.push_back(cm);
				//��cm�е�cover_rs�д������
				ind[it].sol_gene[cm].cover_rs.push_back(i);
				if (t == 1)
				{
					ind[it].sol_gene[cm].flag = 0;   // ��ʾ���в����Ա��Ƴ���ѡ��
					ind[it].sol_gene[cm].is_in_c = 1;
				}
			}
		}
		//���������е�Ȩֵ�Լ�����Ȩ
		for (j = 1; j <= col_num; j++)
		{

			for (i = 1; i <= ind[it].sol_gene[j].cover_rs.size(); i++)
			{
				temp = ind[it].sol_gene[j].cover_rs[i - 1];
				ind[it].sol_gene[j].score += ind[it].sol_meme[temp].weight;
			}
			ind[it].prior[j] = ind[it].sol_gene[j].score * 1.0 / ind[it].sol_gene[j].cost;
		}
		f1.close();
	}
}


/**
  *brief �γɳ�ʼ��
  */
void initialize()
{
	    for (int ii = 0; ii < pop_size; ii++)
	   {
		   int best_array[MAXC];
		   while (!ind[ii].uncover_rs.empty())
		   {
				   srand((unsigned)time(NULL));//���������
			       double max = -DBL_MAX;
			       int cnt = 0;
			       int i, j;
			       for (j = 1; j <= col_num; j++)
				   {
				          if (ind[ii].sol_gene[j].is_in_c)
							  continue;
				          if (max<ind[ii].prior[j])
				          {
					           max = ind[ii].prior[j];
					           best_array[0] = j;
					           cnt = 1;
				          }
				         else if (max == ind[ii].prior[j])
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
				                    change_in(best_array[bestA[i]], &ind[ii]);
			        }
		  }
	  }
}


/**
  *brief ���ܣ�
  *remarks �����ڲ��һ�ûϸ����
  */
void change_in(int j, individual  *Ind)
{
	int temp_r;
	int temp_c;
	bool first;
	//cout<<"MDZZ"<<endl;
	//��Ӧ����λ����Ϊ1
	(*Ind).sol_gene[j].is_in_c = 1;
	//����score
	(*Ind).sol_gene[j].score = -(*Ind).sol_gene[j].score;
	//��������Ȩ
	(*Ind).prior[j] = -(*Ind).prior[j];
	//�����ھӵ�Ȩֵ��״̬, ��ʾ���Ա������ѡ��
	for (int i = 0; i < (*Ind).sol_gene[j].cover_rs.size(); i++)
	{
		//��j�п��Ը��ǵ���
		temp_r = (*Ind).sol_gene[j].cover_rs[i];
		//���ڼ���Ƿ��ǵ�һ�θ���temp_r
		first = true;
		//�ҵ�j�й���temp_r���ھ� �� ����Ȩֵ��״̬
		for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
		{
			//���Ը���temp_r�е���
			temp_c = (*Ind).sol_meme[temp_r].covered_cs[k];
			//�����j��ֱ������
			if (temp_c == j) continue;
			//������ǵ�j�У�����j���ھӸ����ھӵ�����Ȩ
			else
			{
				//���ڴ��ں�ѡ���е���
				if ((*Ind).sol_gene[temp_c].is_in_c == 1)
				{
					first = false;
					(*Ind).sol_gene[temp_c].score += (*Ind).sol_meme[temp_r].weight;
					(*Ind).prior[temp_c] = (*Ind).sol_gene[temp_c].score * 1.0 / (*Ind).sol_gene[temp_c].cost;
				}
				//���ڲ��ں�ѡ����
				else if ((*Ind).sol_gene[temp_c].is_in_c == 0)
				{
					//��j���ǵ�һ�θ��ǵ�ʱ����ھӸ���Ȩֵ
					if (!first)
					{
						(*Ind).sol_gene[temp_c].score -= (*Ind).sol_meme[temp_r].weight;
						(*Ind).prior[temp_c] = (*Ind).sol_gene[temp_c].score * 1.0 / (*Ind).sol_gene[temp_c].cost;
					}
				}
			}
			//����״̬
			(*Ind).sol_gene[temp_c].config = 1;
		}
		//����ǵ�һ�θ���
		if (first)
		{
			for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
			{
				temp_c = (*Ind).sol_meme[temp_r].covered_cs[k];
				if (temp_c == j) continue;
				else
				{
					(*Ind).sol_gene[temp_c].score -= (*Ind).sol_meme[temp_r].weight;
					(*Ind).prior[temp_c] = (*Ind).sol_gene[temp_c].score * 1.0 / (*Ind).sol_gene[temp_c].cost;
				}
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
}

/**
  *brief ���ܣ�
  *remarks �����ڲ��һ�ûϸ����
  */
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

/**
  *brief ������Mind��/�м�Ⱥ�壿���еĽ⣬ִ�б������
  *remarks �����ڲ��һ�ûϸ����
  */
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

/**
  *brief ������Mind��/�м�Ⱥ�壿���еĽ⣬ִ�н�������
  *remarks ������meme����Ƭ�ν��н��档ֱ�ӽ����������indȺ���У�
  */
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


/**
  *brief ��һ�������Ľ⿪ʼ��ִ������ֲ��ѹ���
  *param[in] individual* ind, һ��ָ���Ľ����
  *remarks ����ֲ������Ĳ������ð�����meme����Ƭ���У�����sls()�ĳ�������������õģ���ν��������Ҫָ��������/���ۡ�sls()�ڲ��һ�ûϸ��
  */
void SLS(individual *ind)
{

	int  step = 1;
	int  temp_cost = evaluate(ind);
	while (step <= max_step)
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
  *brief ��������ô���......����cover������ѡ��һ���Ƴ����У�
  *param[in] individual* ind, һ��ָ���Ľ����
  *param[in] int r, ʲô��˼��
  *return maxc, ��ʲô����󸲸ǵ�������/Ԫ������
  *remarks ��������"_out"��ʲô��˼��
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
  *brief ��������ô���......����uncover������ѡ��һ���Ƴ����У�
  *param[in] individual* ind, һ��ָ���Ľ����
  *param[in] int r, ʲô��˼��
  *return maxc, ��ʲô����󸲸ǵ�������/Ԫ������
  *remarks ��������"_in"��ʲô��˼��
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

/**
  *brief �Ƚ������������Ӧ��ֵ
  *remarks ��cmp()�������õĴ����ܶ�ʱ������ʵ�ַ�ʽ�����Ӻܶ�ʱ�俪����Ҫ����ָ��
  */
bool cmp(individual  a, individual b)
{
	return  a.fitness <  b.fitness;
}

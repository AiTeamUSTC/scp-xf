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
  *brief 定义基因片段一每个基因位的结构体，每个基因位代表一个“列”/“子集合”，它能覆盖若干“行”/“元素”
*/
typedef struct gene{
	int config;             ///是否可以加入候选解
	int flag;               ///是否可以移除候选解
	int score;              ///列的分数：值越大表示......，值越小表示......
	int time_stamp;         ///时间戳
	int is_in_c;            ///加入候选解则编码为1，否则为0
	int cost;               ///什么意思？用在哪里？
	vector<int> cover_rs;   /// rows it can cover
}gene;



/**
  *brief 定义基因片段二每个基因位的结构体，一个基因位代表一个“要覆盖”的行（？）
*/
typedef struct meme{
	int weight;             /// 表示一个元素或行被覆盖的“难度”，值越大，表示覆盖越困难
	vector<int> covered_cs; /// cols have can cover this row
}meme;


/**
  *brief 定义种群中每个个体的结构体
*/
typedef struct individual{
	gene sol_gene[MAXC];    ///基因片段1
	meme sol_meme[MAXR];    ///基因片段2
	double prior[MAXC];     ///什么？
	vector<int> uncover_rs; ///未被覆盖的行
	int fitness;            ///适应度的值，怎么计算，表示什么？被选中的列的数目？

}individual;



int row_num;                        ///instance 特征：行数
int col_num;                        ///instance 特征：列数
double Cr = 0.8;                    ///算法参数： 交叉概率？
individual ind[pop_size];           ///算法参数：群体
individual Mind[pop_size];          ///算法参数：？？群体？
vector<int> best_sol;               ///算法执行过程中找到的最优解


///函数的声明
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
  * @brief 主程序
  * @author xufang
  * @param [in] NULL
  * @param [out] NULL
  * @return NULL
  * @note ....todo
  */
int main(int argc, char *argv[])
{
	srand((unsigned)time(NULL));//随机化种子
	char filename[] = "scp41.txt";
	//种群建模
	build_instance(filename);
	//种群初始化
	initialize();
	for (int i = 0; i < pop_size; i++)
	{
		//评估初始种群
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
///主程序结束


/**
  *brief 产生随机数
  */
double random(double low, double high)
{
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}


/**
  *brief 实例建模
  */
void build_instance(char* file)
{
	for (int it = 0; it < pop_size; it++)
	{
		int i, j, k, t, p, cm, temp;
		ifstream f1;
		f1.open(file);
		//读入总的行数以及列数
		f1 >> row_num >> col_num;
		//初始清空uncover_rs
		ind[it].uncover_rs.clear();
		//读入每一列的代价
		for (j = 1; j <= col_num; j++)
		{
			f1 >> k;
			/*如果是要解决单位代价的集合覆盖问题，直接把K赋值为1即可*/
			//ind[it].sol_gene[j].cost = k;
			ind[it].sol_gene[j].cost = 1;
			ind[it].sol_gene[j].config = 1;    //表示都可以加入候选解
			ind[it].sol_gene[j].flag = 1;     // 表示初始状态所有列都可以从候选解移除
			ind[it].sol_gene[j].score = 0;    // 初始化score为0
			ind[it].sol_gene[j].is_in_c = 0;   //初始所有列都不在候选解中
			ind[it].sol_gene[j].time_stamp = 0; // 时间戳
			ind[it].sol_gene[j].cover_rs.clear(); 			//初始清空cover_rs
		}
		// 初始化个体的适应值函数
		ind[it].fitness = 0;
		for (i = 1; i <= row_num; i++)
		{
			//初始清空covered_cs
			ind[it].sol_meme[i].covered_cs.clear();
			//初始化所有的行为未被覆盖
			ind[it].uncover_rs.push_back(i);
		}
		//读入各行列的覆盖关系
		for (i = 1; i <= row_num; i++)
		{
			//初始化每一行的权重
			ind[it].sol_meme[i].weight = 1;
			//读入第i行可以被多少列覆盖
			f1 >> t;
			//读入可以覆盖第i行的列
			for (p = 0; p < t; p++)
			{
				f1 >> cm;
				//在covered_cs中存入cm列
				ind[it].sol_meme[i].covered_cs.push_back(cm);
				//在cm列的cover_rs中存入该行
				ind[it].sol_gene[cm].cover_rs.push_back(i);
				if (t == 1)
				{
					ind[it].sol_gene[cm].flag = 0;   // 表示该列不可以被移除候选解
					ind[it].sol_gene[cm].is_in_c = 1;
				}
			}
		}
		//计算所有列的权值以及优先权
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
  *brief 形成初始解
  */
void initialize()
{
	    for (int ii = 0; ii < pop_size; ii++)
	   {
		   int best_array[MAXC];
		   while (!ind[ii].uncover_rs.empty())
		   {
				   srand((unsigned)time(NULL));//随机化种子
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
  *brief 功能？
  *remarks 函数内部我还没细看。
  */
void change_in(int j, individual  *Ind)
{
	int temp_r;
	int temp_c;
	bool first;
	//cout<<"MDZZ"<<endl;
	//相应基因位编码为1
	(*Ind).sol_gene[j].is_in_c = 1;
	//更新score
	(*Ind).sol_gene[j].score = -(*Ind).sol_gene[j].score;
	//更新优先权
	(*Ind).prior[j] = -(*Ind).prior[j];
	//更新邻居的权值和状态, 表示可以被加入候选解
	for (int i = 0; i < (*Ind).sol_gene[j].cover_rs.size(); i++)
	{
		//第j列可以覆盖的行
		temp_r = (*Ind).sol_gene[j].cover_rs[i];
		//用于检测是否是第一次覆盖temp_r
		first = true;
		//找第j列关于temp_r的邻居 ， 更新权值和状态
		for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
		{
			//可以覆盖temp_r行的列
			temp_c = (*Ind).sol_meme[temp_r].covered_cs[k];
			//如果是j列直接跳过
			if (temp_c == j) continue;
			//如果不是第j列，就是j的邻居更新邻居的优先权
			else
			{
				//对于处于候选解中的列
				if ((*Ind).sol_gene[temp_c].is_in_c == 1)
				{
					first = false;
					(*Ind).sol_gene[temp_c].score += (*Ind).sol_meme[temp_r].weight;
					(*Ind).prior[temp_c] = (*Ind).sol_gene[temp_c].score * 1.0 / (*Ind).sol_gene[temp_c].cost;
				}
				//对于不在候选解内
				else if ((*Ind).sol_gene[temp_c].is_in_c == 0)
				{
					//且j不是第一次覆盖的时候的邻居更新权值
					if (!first)
					{
						(*Ind).sol_gene[temp_c].score -= (*Ind).sol_meme[temp_r].weight;
						(*Ind).prior[temp_c] = (*Ind).sol_gene[temp_c].score * 1.0 / (*Ind).sol_gene[temp_c].cost;
					}
				}
			}
			//更新状态
			(*Ind).sol_gene[temp_c].config = 1;
		}
		//如果是第一次覆盖
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
			//将temp_r从uncover_rs中删除
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
  *brief 功能？
  *remarks 函数内部我还没细看。
  */
void change_out(int j, individual * ind)
{
	int temp_r;
	int temp_c;
	bool first;
	//相应基因位编码为0
	(*ind).sol_gene[j].is_in_c = 0;
	//更新score
	(*ind).sol_gene[j].score = -(*ind).sol_gene[j].score;
	//更新优先权
	(*ind).prior[j] = -(*ind).prior[j];
	(*ind).sol_gene[j].config = 0;
	//更新邻居的状态, 表示可以被加入候选解
	for (int i = 0; i < (*ind).sol_gene[j].cover_rs.size(); i++)
	{
		//第j列可以覆盖的行
		temp_r = (*ind).sol_gene[j].cover_rs[i];
		first = true;
		//找第j列关于temp_r的邻居
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
		//当候选解中依然存在可以覆盖temp_r的列时
		if (!first)
		{
			//找出可以覆盖temp_r的所有列
			for (int k = 0; k < (*ind).sol_meme[temp_r].covered_cs.size(); k++)
			{
				temp_c = (*ind).sol_meme[temp_r].covered_cs[k];
				if (temp_c == j) continue;
				//更新候选解外的邻居的权值
				else if (temp_c != j && (*ind).sol_gene[temp_c].is_in_c == 0)
				{
					(*ind).sol_gene[temp_c].score -= (*ind).sol_meme[temp_r].weight;
					(*ind).prior[temp_c] = (*ind).sol_gene[temp_c].score * 1.0 / (*ind).sol_gene[temp_c].cost;
				}
				//更新候选解内的邻居的权值
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
			//更新邻居的权值
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
  *brief 对所有Mind（/中间群体？）中的解，执行变异操作
  *remarks 函数内部我还没细看。
  */
void mutation()
{

	//通过突变的操作生成新的种群
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
	//通过突变操作生成新的种群
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
  *brief 对所有Mind（/中间群体？）中的解，执行交叉算子
  *remarks 仅仅对meme基因片段进行交叉。直接将结果保存在ind群体中？
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
  *brief 从一个给定的解开始，执行随机局部搜过程
  *param[in] individual* ind, 一个指定的解个体
  *remarks 随机局部搜索的参数设置包含在meme基因片段中，但是sls()的超参数是如何设置的？所谓超参数主要指搜索步数/代价。sls()内部我还没细看
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
			continue;//这个地方不用while是为了保证初始阶段change_out的列的time_stamp不同
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
		//未被覆盖的行的权值增加，相应的可以覆盖该行的列的score和优先权增加
		for (int i = 0; i < (*ind).uncover_rs.size(); i++)
		{
		int ur = (*ind).uncover_rs[i];
		(*ind).sol_meme[ur].weight += 1;
		for (int j = 0; j < (*ind).sol_meme[ur].covered_cs.size(); j++)
		{
		//可以覆盖该行的列
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
  *brief 计算解/个体的适应度值
  *param[in] individual* ind, 一个指定的解个体
  *return int fitness, 被选中的列的数目
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
  *brief 这个函数用处是......？从cover集合中选择一个移除的列？
  *param[in] individual* ind, 一个指定的解个体
  *param[in] int r, 什么意思？
  *return maxc, 是什么？最大覆盖掉的行数/元素数？
  *remarks 函数名中"_out"是什么意思？
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
  *brief 这个函数用处是......？从uncover集合中选择一个移除的列？
  *param[in] individual* ind, 一个指定的解个体
  *param[in] int r, 什么意思？
  *return maxc, 是什么？最大覆盖掉的行数/元素数？
  *remarks 函数名中"_in"是什么意思？
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
  *brief 比较两个个体的适应度值
  *remarks 当cmp()函数调用的次数很多时，这种实现方式会增加很多时间开销，要改用指针
  */
bool cmp(individual  a, individual b)
{
	return  a.fitness <  b.fitness;
}

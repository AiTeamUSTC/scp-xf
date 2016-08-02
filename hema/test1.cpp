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

/********定义基因片段一每个基因位的结构体，每个基因位包含很多信息，可以覆盖的行，该列的代价，时间戳等等********/
typedef struct gene{
	int config;  //是否可以加入候选解 
	int flag;    //是否可以移除候选解 
	int score; //列的分数，初始值为它所能覆盖的行的个数，在单位代价集合覆盖问题中，score的值越大的列，越优先被选入候选解
	int time_stamp; //时间戳 
	int is_in_c;    // 加入候选解则编码为1，否则为0 
	int cost;  //每一列的代价，实例的中前2000个数据，表示的就是每一行的代价，单位代价集合覆盖问题中，值统一为1
	vector<int> cover_rs;// 存放该列可以覆盖的行
}gene;


/********定义基因片段二每个基因位的结构体，每个基因位代表一行，并包含行的权值属性以及每一行可以被哪些列覆盖********/
typedef struct meme{
	int weight;  // 每一行的权重，在求解过程中不断变化，越难覆盖到的行，权值越大，从而增加其被覆盖的可能性
	vector<int> covered_cs; // 存放可以覆盖该行的列
}meme;


/********定义种群中每个个体的结构体********/
typedef struct individual{
	gene sol_gene[MAXC];  //基因片段一
	meme sol_meme[MAXR]; //基因片段二
	double prior[MAXC];       //用来定义基因片段一中相应基因位被选入候选解的优先权,单位代价集合覆盖问题中，与score保持一致
	vector<int> uncover_rs;   //存放每个解没有覆盖到的行
	int fitness;                       // 简单表示候选解的大小
}individual;


/********定义全局变量************/
int row_num;//行数
int col_num;//列数
double Cr = 0.8; // 交叉率
individual ind[pop_size];//原种群
individual Mind[pop_size];//突变种群
/*********函数的定义***********************/
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


/*********主程序*******************/
int main(int argc, char *argv[])
{
	 //srand((unsigned)time(NULL));//随机化种子 
	char filename[] = "scp41.txt";
	//种群建模 
	build_instance(filename);
	srand(40);
	for (int i = 0; i < pop_size; i++)
	{
		//种群初始化 
		  initialize(&ind[i]);
		//评估初始种群 
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
/*********主程序结束*******************/


/****产生随机数 ***/
double random(double low, double high)
{
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}

/******实例建模********/
void build_instance(char* file)
{
	for (int it = 0; it < pop_size; it++)
	{
		int i, j, k, t, p, cm, temp;
		ifstream f1;
		f1.open(file);
		f1 >> row_num >> col_num;		//读入总的行数以及列数 
		
		/**************初始化阶段*************************/ 
		ind[it].uncover_rs.clear();		//初始清空uncover_rs 
		for (i = 1; i <= row_num; i++)
		{
			ind[it].sol_meme[i].covered_cs.clear();	         //初始清空所有行的covered_cs
		}
		for (j = 1; j <= col_num; j++)
		{
			f1 >> k;		//读入每一列的代价 
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
		ind[it].fitness = 0; 		// 初始化个体的适应值函数 
		
	   /****************开始读入相应的行列关系********************/ 
		for (i = 1; i <= row_num; i++)
		{
			ind[it].uncover_rs.push_back(i);			//初始化所有的行为未被覆盖
			ind[it].sol_meme[i].weight = 1;			//初始化每一行的权重 
			f1 >> t;			//读入第i行可以被多少列覆盖，这里是t列  
			if (t == 1)// 如果存在某一行只能被一列覆盖，则将该列固定在候选解中 
			{
				ind[it].sol_gene[cm].flag = 0;   // 表示该列不可以被移出候选解
				ind[it].sol_gene[cm].is_in_c = 1; // 表示该列在候选解中 
				ind[it].uncover_rs.pop_back(); //相应被覆盖的那一行移出uncover_rs
			}
			for (p = 0; p < t; p++)
			{
				f1 >> cm;			//读入可以覆盖第i行的列
				ind[it].sol_meme[i].covered_cs.push_back(cm);				//在covered_cs中存入cm列 
				ind[it].sol_gene[cm].cover_rs.push_back(i);				//在cm列的cover_rs中存入该行 
			}
		}
		//计算所有列的权值以及优先权 
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

/*****形成初始解 ******/
void initialize(individual * ind)
{
		   int best_array[MAXC]; 
		   double max;
		   int cnt; 
		   int i, j;
		   while (!(*ind).uncover_rs.empty())
		   {
				   //srand((unsigned)time(NULL));//随机化种子

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
/*change_in函数，主要对基因片段一中基因位的is_in_c属性做改变，is_in_c从0到1表示将相应列加入候选解，函数中也包括对邻居的优先权更新*/ 
void change_in(int j, individual  *Ind)
{
	   int temp_r, temp_c, btemp_c, temp, count;
 	   (*Ind).sol_gene[j].is_in_c = 1;	   //相应基因位属性编码为1,表示将j列加入候选解
	   (*Ind).sol_gene[j].score = -(*Ind).sol_gene[j].score; 	   // score置为相反数  
	   (*Ind).prior[j] = -(*Ind).prior[j]; 	   //相应的更新优先权 
	   	/*找出第j列的邻居，并更新邻居的优先权以及状态*/ 
       for (int i = 0; i < (*Ind).sol_gene[j].cover_rs.size(); i++)
	   {
		   temp_r = (*Ind).sol_gene[j].cover_rs[i];     //第j列可以覆盖的行 
		   count = 0;  			//用于检测是否是第一次覆盖temp_r 
		   /*找第j列关于temp_r的邻居 ， 更新权值和状态 */
		   for (int k = 0; k < (*Ind).sol_meme[temp_r].covered_cs.size(); k++)
		   {
			   temp_c = (*Ind).sol_meme[temp_r].covered_cs[k]; 				 //可以覆盖temp_r行的列 
			   if (temp_c == j)   continue;				 //如果是j列直接跳过 
			   /*判断temp_c列是否在候选解中*/
			   if ((*Ind).sol_gene[temp_c].is_in_c == 1)
			   {
				   btemp_c = temp_c; //记录下在候选解中的邻居 
				   count++;
			   }

			   (*Ind).sol_gene[temp_c].config = 1;  		//更新所有邻居状态
		   }
		   /*count为0表示是temp_r第一次覆盖，相应的第j列关于temp_r的所有邻居的score需减去temp_r行的权值，并更新优先权*/ 
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
		   /*在候选解中的邻居score项应该是负值，所以每多一个新邻居进入候选解，相应的temp_r的权值的影响就需要被抹去*/
		   /*至于为什么只抹去第一个，论文解释是当temp_r被多次覆盖时，temp_r的权值对列score的影响可以不计*/ 
		   else if (count == 1)
		   {
			   (*Ind).sol_gene[btemp_c].score += (*Ind).sol_meme[temp_r].weight;
			   (*Ind).prior[btemp_c] = (*Ind).sol_gene[btemp_c].score * 1.0 / (*Ind).sol_gene[btemp_c].cost;
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
/*change_out函数，主要对基因片段一中基因位的is_in_c属性做改变，is_in_c从1到0表示将相应列移出候选解，函数中也包括对邻居的优先权更新*/ 
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
/*突变操作，主要是形成突变种群，突变发生在基因片段二，方便后面利用突变种群和原种群做重组主要参考了一般DE算法的突变重组操作*/ 
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
/*重组操作，随机对ind中的解的某些列，执行交叉算子，仅仅对meme基因片段进行交叉，直接将结果保存在ind群体中*/
/*这里突变重组操作仅仅是起到一个种群间信息交流的作用，选择以及对解的改进都集中在SLS里了*/ 
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
 /*从一个给定的解开始，执行随机局部搜过程
  * individual* ind, 一个指定的解个体
  *随机局部搜索的参数设置包含在meme基因片段中，但是sls()的超参数在实验调试阶段为缩短时间，只是取1000；
  *进一步改进解的质量时会设置的更大一些 ，RWLS是30000000 
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
  *brief 从候选解中选择score最大的列移出 
  *param[in] individual* ind, 一个指定的解个体
  *return maxc，候选解中被选择移出的列
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
  *brief 从非候选解的列中选择score最大的列移入候选解 
  *param[in] individual* ind, 一个指定的解个体
  *return maxc，被选择移入候选解的列
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
/*写这个函数目的是想种群排个序，选择最优的个体参与到突变重组中，目前还没用到*/ 
bool cmp(individual  *a, individual *b)
{
	return  (*a).fitness <  (*b).fitness;
}



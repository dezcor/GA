#ifndef _GA_H_
#define _GA_H_
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <exception>

using namespace std;

namespace GA
{
	class Individuo
	{
	private:
		vector<unsigned int> Cromosoma;
		float Fitness;
		float Objetivo;
		float ProbSel;
		vector<float> valor;

	public:
		//___________Get____________
		unsigned int Getsize();
		float GetFitness();
		float GetObjetivo();
		float GetValor(const int i);
		float GetProbSel();
		vector<unsigned int>& GetCromosoma();
		//___________Set____________
		unsigned int &operator[](const int Index);
		void SetNumeroGenes(const unsigned int size);
		void SetNumeroBitGen(const vector<unsigned int> &NumBitGen);
		void setValor(const unsigned int Index, const float valor);
		void SetObjetivo(const float obj);
		void SetFitness(float fitness);
		void SetProbSel(float pbSel);
		void SetCromosoma(vector<unsigned int>& Cromosoma);

		Individuo();
		~Individuo() {}
	};

	class Poblacion
	{
	protected:
		unsigned int NumeroDeGenes;
		vector<unsigned int> NumeroBitGet;
		vector<Individuo> individuos;
		vector<Individuo> NewIndividuos;
		vector<unsigned int> idSelect;
		unsigned int SizeCromosoma;
		unsigned int SizePoblacion;
		float ObjetivoMinimo;
		float SumFitness;
		float Pc;
		float Pm;
		int idMejor;
		int MaxError;
		bool isNumberReal;
		std::vector<float> Limites;
	public:
		float Error;
		//____________Get____________
		unsigned int GetNumGenes();
		unsigned int GetNumBitGet(const unsigned int Index);
		unsigned int GetsizeCromosoma();
		unsigned int GetsizePoblacion();
		int GetidMejor();
		Individuo &operator[](const unsigned int Index);
		bool isNumberR()
		{
			return isNumberReal;
		}
		//____________Set____________
		void SetNumGenes(const unsigned int NG);
		void SetNumBitGen(const vector<unsigned int> &NumeroBitGet);
		void resizePoblacion(const unsigned int NumeroDeGenes);
		void setProbabilidadCruza(float pc);
		void setProbabilidadMuta(float pm);
		void SetisNumberReal(bool is)
		{
			isNumberReal= is;
		}
		void SetLimites(std::vector<float> const & Lim)
		{
			Limites = Lim;
			isNumberReal = true;
		}
		Poblacion();
		virtual ~Poblacion() {}

		//_____________________________Algoritmo
		void Inicializa();
		void Ruleta();
		void Cruza();
		void Muta();
		void Evaluacion();
		void Fitness();
		virtual void FuncionObjetivo(Individuo& individuo);
		void ConservarMejor();
		void ActualizarPoblacion();
		void MostrarIndividuo(const unsigned int Index);
		void MostrarPoblacion();
		//********************Decodificar***********************//
		void DecodificarToInt(Individuo& individuo);
		void Run(int64_t Iteraciones, float Error);
		//___________________________End
	};
}; // namespace GA

#endif //_GA_H_
#include "GA.h"

namespace GA
{
	/***********************************
	* Class Individuo
	*
	*
	***********************************/

	Individuo::Individuo()
	{
		Fitness = 0;
		this->Objetivo = 0;
		this->ProbSel = 0;
	}
	//___________Get____________
	float Individuo::GetFitness()
	{
		return Fitness;
	}

	float Individuo::GetObjetivo()
	{
		return Objetivo;
	}

	unsigned int Individuo::Getsize()
	{
		return (unsigned int)Cromosoma.size();
	}

	float Individuo::GetValor(int Index)
	{
		if (Index >= valor.size())
			return 0;

		return valor[Index];
	}
	float Individuo::GetProbSel()
	{
		return this->ProbSel;
	}

	vector<unsigned int>& Individuo::GetCromosoma()
	{
		return Cromosoma;
	}
	//___________Set____________
	unsigned int &Individuo::operator[](const int Index)
	{
		return Cromosoma[Index];
	}

	void Individuo::SetNumeroGenes(const unsigned int size)
	{
		valor.resize(size);
	}

	void Individuo::SetNumeroBitGen(const vector<unsigned int> &NumBitGen)
	{
		SetNumeroGenes(NumBitGen.size());
		int sum = 0;
		for (unsigned int Get : NumBitGen)
			sum += Get;
		Cromosoma.resize(sum);
	}

	void Individuo::setValor(const unsigned int Index, const int Valor)
	{
		this->valor[Index] = Valor;
	}

	void Individuo::SetObjetivo(const float obj)
	{
		this->Objetivo = obj;
	}

	void Individuo::SetFitness(float fitness)
	{
		this->Fitness = fitness;
	}

	void Individuo::SetProbSel(float pbSel)
	{
		this->ProbSel = pbSel;
	}

	void Individuo::SetCromosoma(vector<unsigned int>& Cromosoma)
	{
		this->Cromosoma = Cromosoma;
	}

	/******************************************
	* Class Poblacion
	******************************************/
	Poblacion::Poblacion()
	{
		NumeroDeGenes = 0;
		SizeCromosoma = 0;
		SizePoblacion = 0;
		Pm = 0;
		Pc = 0;
		MaxError = 50;
		idMejor = 0;
		srand(time(NULL));
	}

	unsigned int Poblacion::GetNumGenes()
	{
		return NumeroDeGenes;
	}

	unsigned int Poblacion::GetNumBitGet(const unsigned int Index)
	{
		if (Index >= NumeroBitGet.size())
			return 0;
		return NumeroBitGet[Index];
	}

	unsigned int Poblacion::GetsizeCromosoma()
	{
		return SizeCromosoma;
	}

	unsigned int Poblacion::GetsizePoblacion()
	{
		return individuos.size();
	}

	Individuo &Poblacion::operator[](const unsigned int Index)
	{
		return individuos[Index];
	}

	int Poblacion::GetidMejor()
	{
		return idMejor;
	}

	void Poblacion::SetNumGenes(const unsigned int NG)
	{
		NumeroDeGenes = NG;
	}

	void Poblacion::SetNumBitGen(const vector<unsigned int> &NumeroBitGet)
	{
		this->NumeroBitGet = NumeroBitGet;
		int sum = 0;
		for (int i = 0; i < NumeroBitGet.size(); i++)
			sum += NumeroBitGet[i];
		SizeCromosoma = sum;
		NumeroDeGenes = NumeroBitGet.size();
	}

	void Poblacion::resizePoblacion(const unsigned int size)
	{
		individuos.resize(size);
		NewIndividuos.resize(size);
		SizePoblacion = size;
		idSelect.resize(size);
	}

	void Poblacion::setProbabilidadCruza(float pc)
	{
		Pc = pc;
	}

	void Poblacion::setProbabilidadMuta(float pm)
	{
		Pm = pm;
	}

	//_____________________________Algoritmo Virtuales
	void Poblacion::Inicializa()
	{
		unsigned int i, j;
		for (i = 0; i < SizePoblacion; i++)
		{
			individuos[i].SetNumeroBitGen(this->NumeroBitGet);
			NewIndividuos[i].SetNumeroBitGen(this->NumeroBitGet);
			for (j = 0; j < SizeCromosoma; j++)
			{
				individuos[i][j] = rand() % 2;
			}
		}
		Error = 100;
	}

	void Poblacion::Ruleta()
	{
		/*******************Insertar Codigo*************/
		float maxR = 1.0 / SumFitness;
		for (int i = 0; i < this->GetsizePoblacion(); i++)
		{
			individuos[i].SetProbSel(individuos[i].GetFitness() * maxR);
		}
		float sum, aux;
		int i, j;
		for (i = 0, j = 0, sum = 0; i < this->GetsizePoblacion(); i++, j = 0, sum = 0)
		{
			aux = (double)rand() / RAND_MAX;
			do
			{
				sum += individuos[j].GetProbSel();
				j++;
			} while (aux > sum && j < GetsizePoblacion());
			idSelect[i] = j - 1;
		}
	}

	void Poblacion::Evaluacion()
	{
		/*******************Insertar Codigo*************/
		float miniObj = INT64_MAX;
		for (int i = 0; i< this->GetsizePoblacion(); i++)
		{
			FuncionObjetivo(individuos[i]);
			if (individuos[i].GetObjetivo() < miniObj)
			{
				miniObj = individuos[i].GetObjetivo();
			}
		}
		this->ObjetivoMinimo = miniObj;
		Fitness();
	}

	void Poblacion::Cruza()
	{
		/*******************Insertar Codigo*************/
		int RangoCruza; //= (int)((double)rand() / RAND_MAX * (SizeCromosoma - 2) + 1);
		float random = 0;
		for (int i = 0; i < SizePoblacion - 1; i += 2)
		{
			random = (double)rand() / RAND_MAX;
			if (random < Pc)
			{
				RangoCruza = (int)((double)rand() / RAND_MAX * (SizeCromosoma - 2) + 1);
				for (int j = 0; j < RangoCruza; j++)
				{
					NewIndividuos[i][j] = individuos[idSelect[i]][j];
					NewIndividuos[i + 1][j] = individuos[idSelect[i + 1]][j];
				}
				for (int j = RangoCruza; j < SizeCromosoma; j++)
				{
					NewIndividuos[i][j] = individuos[idSelect[i + 1]][j];
					NewIndividuos[i + 1][j] = individuos[idSelect[i]][j];
				}
			}
			else
			{
				NewIndividuos[i].SetCromosoma(individuos[i].GetCromosoma());
				NewIndividuos[i + 1].SetCromosoma(individuos[i + 1].GetCromosoma());
			}
		}
	}

	void Poblacion::Fitness()
	{
		SumFitness = 0;
		for (int i = 0; i < this->GetsizePoblacion(); i++)
		{
			individuos[i].SetFitness(individuos[i].GetObjetivo() - this->ObjetivoMinimo);
			SumFitness += individuos[i].GetFitness();
		}
		/*******************Insertar Codigo*************/
	}

	void Poblacion::FuncionObjetivo(Individuo &individuo)
	{
		/*******************Insertar Codigo*************/
		DecodificarToInt(individuo);
		int x = individuo.GetValor(0);
		int y = individuo.GetValor(1);
		individuo.SetObjetivo(50 - (x - 5) * (x - 5) - (y - 5) * (y - 5));
	}

	void Poblacion::Run(int64_t MaxGeneraciones, float Error)
	{
		/*******************Insertar Codigo*************/
		int generacion = 1;
		Inicializa();
		Evaluacion();
		MostrarPoblacion();
		while (generacion < MaxGeneraciones && this->Error> Error)
		{
			Ruleta();
			Cruza();
			Muta();
			ConservarMejor();
			ActualizarPoblacion();
			Evaluacion();
			printf("Generacion %d\n\n", generacion);
			//MostrarPoblacion();
			generacion++;
		}

	}

	void Poblacion::Muta()
	{
		float random;
		for (int i = 0; i < SizePoblacion; i++)
		{
			for (int j = 0; j < SizeCromosoma; j++)
			{
				random = (double)rand() / RAND_MAX;
				if (random < Pm)
				{
					NewIndividuos[i][j] = !NewIndividuos[i][j];
				}
			}
		}
		/*******************Insertar Codigo*************/
	}
	void Poblacion::ConservarMejor()
	{
		for (int i = 0; i < SizePoblacion; i++)
		{
			if (individuos[this->idMejor].GetObjetivo() < individuos[i].GetObjetivo())
			{
				idMejor = i;
			}
		}
		Error = fabs(MaxError - individuos[this->idMejor].GetObjetivo());
	}

	void Poblacion::ActualizarPoblacion()
	{
		for (int i = 0; i < SizePoblacion; i++)
		{
			if (i != idMejor)
			{
				individuos[i].SetCromosoma(NewIndividuos[i].GetCromosoma());
			}

		}
	}

	void Poblacion::DecodificarToInt(Individuo &individuo)
	{
		int Init = 0, next = 0, Numero;
		for (int i = 0; i < this->NumeroDeGenes; i++)
		{
			next += NumeroBitGet[i];
			Numero = 0;
			for (int j = Init, k = 0; j < next; j++, k++)
			{
				Numero = Numero | ((individuo[j] == 1) ? 1 << k : 0);
			}
			Init = next;
			individuo.setValor(i, Numero);
		}
	}

	void Poblacion::MostrarIndividuo(const unsigned int Index)
	{
		for (int j = SizeCromosoma - 1; j >= 0; j--)
		{
			printf("%d", individuos[Index][j]);
		}
		printf("\n");
		for (int j = NumeroDeGenes - 1; j >= 0; j--)
		{
			printf("%02.02f ", individuos[Index].GetValor(j));
		}
		printf("\t%f\n", individuos[Index].GetObjetivo());
	}
	void Poblacion::MostrarPoblacion()
	{
		for (int i = 0; i < SizePoblacion; i++)
		{
			printf("Individuo: %d\n", i);
			MostrarIndividuo(i);
		}
	}
	//___________________________End

}; // namespace GA
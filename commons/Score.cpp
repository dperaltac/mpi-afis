#include <string>

#include "Score.h"
using namespace std;

Score::Score() : score(-1), relindex(-1), absindex(-1), process(-1)
{
	strcpy(id, "");
}

Score::Score(const Score &s) :
	score(s.score),
	relindex(s.relindex),
	absindex(s.absindex),
	process(s.process)
{
	strcpy(id, s.id);
}

Score::Score(float s, int reli, int absi, int proc, const std::string &id) :
	score(s),
	relindex(reli),
	absindex(absi),
	process(proc)
{
	strcpy(this->id, id.c_str());
}

Score & Score::operator=(const Score &s)
{
	if(&s != this)
	{
		score = s.score;
		relindex = s.relindex;
		absindex = s.absindex;
		process = s.process;
		strcpy(id, s.id);
	}

	return *this;
}


#ifdef MPI_VERSION
#include "mpi.h"

MPI::Datatype Score::getDatatype()
{
	MPI::Datatype obj_mpi;
	Score testscore;

	MPI::Datatype types[5] = {MPI::FLOAT, MPI::INT, MPI::INT, MPI::INT, MPI::CHAR};
	int block_sizes[5] = {1, 1, 1, 1, IDLENGTH};
	MPI::Aint offset[5];

	offset[0] = MPI::Get_address(&testscore.score)    - MPI::Get_address(&testscore);
	offset[1] = MPI::Get_address(&testscore.relindex) - MPI::Get_address(&testscore);
	offset[2] = MPI::Get_address(&testscore.absindex) - MPI::Get_address(&testscore);
	offset[3] = MPI::Get_address(&testscore.process)  - MPI::Get_address(&testscore);
	offset[4] = MPI::Get_address(testscore.id)        - MPI::Get_address(&testscore);

	obj_mpi = MPI::Datatype::Create_struct(5, block_sizes, offset, types);
	obj_mpi.Commit();

	return obj_mpi;
}

#endif

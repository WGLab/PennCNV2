#include<CL/cl.h>
using namespace std;

const char * clError (cl_int rc);
void clSafe (cl_int rc, string functionname);

inline void checkErr(cl_int err, const char * name)
{
    if (err != CL_SUCCESS) {
        std::cerr << "ERROR: " << name
                 << " (" << err << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
}


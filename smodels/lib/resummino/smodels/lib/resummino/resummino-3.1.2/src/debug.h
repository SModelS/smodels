#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#ifndef NDEBUG

#include <stdio.h>
#include <string>
#include <assert.h>
#include <string.h>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>


using namespace std;

const static char* debugfile = "debug.csv";
const static char* debugtablefile = "debug_table.csv";
const int __debug_rm_return__ = remove(debugfile);
const int __debug_rm_table_return__ = remove(debugtablefile);
extern map<string,bool> lines;
extern map<string,string> cur_table_line;


#define __FILENAME__ (strrchr("/" __FILE__, '/')+1)
// "\\"for windows
#define _DEBUG_FILE(key,value) \
{ \
	std::ofstream ost {debugfile, std::ios_base::app}; \
	ost << key << ";" << setprecision(15) << value << std::endl; \
	ost.close(); \
} 
#define _DEBUG_MSG(format, args...)                                  \
    {                                                                   \
            printf("%s:%i: " format "\n", __FILENAME__,__LINE__ , ##args);\
            std::cout << std::flush;                                    \
    }
#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)
#define _DEBUG_MSG1(format, args...)       \
if(lines.find("DEBUG_MSG1" __FILE__ S__LINE__)==lines.end()){_DEBUG_MSG(format,##args);lines["DEBUG_MSG1" __FILE__ S__LINE__]=true;}
#define _DEBUG_FILE1(key,value)       \
if(lines.find("DEBUG_FILE1" __FILE__ S__LINE__)==lines.end()){_DEBUG_FILE(key,value);lines["DEBUG_FILE1" __FILE__ S__LINE__]=true;}
#define _DEBUG_ASSERT(expr) { assert(expr); }

#define _DEBUG_TABLE_FLUSH_HEADER()  \
{ \
	std::ofstream ost {debugtablefile, std::ios_base::app}; \
	for (std::map<string, string>::iterator i = cur_table_line.begin(); i != cur_table_line.end(); i++) \
	{\
    		ost << i->first << ";"; \
	}\
	ost << std::endl; \
	ost.close(); \
}

#define _DEBUG_TABLE_FLUSH_LINE()  \
{ \
	std::ofstream ost {debugtablefile, std::ios_base::app}; \
	for (std::map<string, string>::iterator i = cur_table_line.begin(); i != cur_table_line.end(); i++) \
	{\
    		ost  << i->second << ";"; \
	}\
	ost << std::endl; \
	ost.close(); \
	cur_table_line.clear();\
}


#define _DEBUG_TABLE(key,value)  \
{\
	std::stringstream ss; \
	ss << scientific <<setprecision(30)<< value; \
	cur_table_line[key]=ss.str();\
}

#else
#define _DEBUG_MSG(format, args...)
#define _DEBUG_MSG1(format, args...)
#define _DEBUG_FILE(format, args...)
#define _DEBUG_FILE1(format, args...)
#define _DEBUG_TABLE(format, args...)
#define _DEBUG_TABLE_FLUSH_LINE(format, args...)
#define _DEBUG_TABLE_FLUSH_HEADER(format, args...)
#define _DEBUG_ASSERT(expr)
#endif
#endif

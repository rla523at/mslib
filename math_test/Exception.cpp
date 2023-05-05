#include "gtest\gtest.h"
#include "../math/src/Exception.h"

#ifdef _DEBUG

TEST(Exception, File_Name)
{
	EXPECT_STREQ(FILE_NAME.data(), "Exception.cpp");
}
//TEST(Exception, Exception)
//{
//	const char* message = "error message for testing";
//	int exception_line_num = 0;
//
//	try
//	{
//		exception_line_num = __LINE__ + 1;
//		EXCEPTION(message);
//	}
//	catch(const std::exception& exception)
//	{
//		std::ostringstream os;
//		os << "\n==============================EXCEPTION========================================\n";
//		os << "File\t\t: " << FILE_NAME << "\n";
//		os << "Function\t: " << __FUNCTION__ << "\n";
//		os << "Line\t\t: " << exception_line_num << "\n";
//		os << "Message\t\t: " << message << "\n";
//		os << "==============================EXCEPTION========================================\n\n";
//
//		EXPECT_EQ(exception.what(), os.str());		
//	}
//}

#else
#endif
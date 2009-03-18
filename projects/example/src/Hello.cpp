#include "Hello.hpp"
#include "Exception.hpp"

Hello::Hello(const std::string& rMessage)
    : mMessage(rMessage)
{
}

std::string Hello::GetMessage()
{
    return mMessage;
}

void Hello::Complain(const std::string& rComplaint)
{
    EXCEPTION(rComplaint);
}

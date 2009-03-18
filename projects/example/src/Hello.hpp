#ifndef HELLO_HPP_
#define HELLO_HPP_

#include <string>

class Hello
{
private:
    std::string mMessage;

public:
    Hello(const std::string& rMessage);
    
    std::string GetMessage();
    
    void Complain(const std::string& rComplaint);
};

#endif /*HELLO_HPP_*/

/* 
 * File:   myExceptions.h
 * Author: tg
 *
 * Created on November 8, 2013, 10:16 AM
 */

#ifndef MYEXCEPTIONS_H
#define MYEXCEPTIONS_H

#include <sstream>

namespace myExcepts {

    class myException {
    protected:
        const std::string msg_;
    public:

        myException(const std::string& text) : msg_(text) {
        }

        ~myException() throw () {
        }
        virtual const std::string what() const throw () = 0;
        std::string msg() const { return msg_; }
    };

    class HeaderEntry : public myException {
    public:

        HeaderEntry(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error: SFRM Header Entry " << msg_
                    << " could not be interpreted\n";
            return msg.str();
        }
    };

    class HeaderFmt : public myException {
    public:

        HeaderFmt(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error: SFRM Invalid Header Format ";
            return msg.str();
        }
    };

    class Usage : public myException {
    public:

        Usage(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error: Invalid command line: " << msg_ << "\n";
            return msg.str();
        }
    };

    class Format : public myException {
    public:

        Format(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error:  File Format Error ";
            return msg.str();
        }
    };

    class FileIO : public myException {
    public:

        FileIO(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error:  Cannot open file " << msg_;
            return msg.str();
        }
    };

    /**
     * \class FileMissing
     * \brief Thrown when expected file could not be opened
     */
    class FileMissing : public myException {
    public:

        FileMissing(const std::string& text) : myException(text) {
        }

        const std::string what() const throw () {
            std::ostringstream msg;
            msg << " *** Error: File " << msg_
                    << " cannot be opened\n";
            return msg.str();
        }
    };
    
    /**
     * \class Programming
     * \brief thrown in case of poor programming leads to undesired call
     */
    class Programming: public myException {
    public:
        Programming(const std::string& text): myException(text) {}
        const std::string what() const throw() {
            std::ostringstream msg;
            msg << " *** Error: Programming leads do undesired behaviour"
                    << msg_ << '\n';
            return msg.str();
        }
    };
    namespace FCF {

        class ListCode : public myException {
        public:

            ListCode(const std::string& text) : myException(text) {
            }

            const std::string what() const throw () {
                std::ostringstream msg;
                msg << " *** Error: fcf-file not LIST code 6\n"
                        << "       " << msg_ << "\n";
                return msg.str();
            }
        };

        class Symop : public myException {
        public:

            Symop(const std::string& text) : myException(text) {
            }

            const std::string what() const throw () {
                std::ostringstream msg;
                msg << " *** Error: The string \""
                        << msg_
                        << "\" cannot be interpreted as Symmetry operator\n";
                return msg.str();
            }
        };

    }

}

#endif /* MYEXCEPTIONS_H */


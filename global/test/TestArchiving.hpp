/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTARCHIVING_HPP_
#define TESTARCHIVING_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <cxxtest/TestSuite.h>

#include "OutputFileHandler.hpp"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/list.hpp>

#include <boost/serialization/export.hpp>

// see http://www.boost.org/libs/serialization/doc/index.html
class ParentClass;

class ChildClass
{
public:
    unsigned mTag;
    ParentClass* mpParent;
    ChildClass() : mTag(1)
    {
    }
    void SetParent(ParentClass* pParent)
    {
        mpParent = pParent;
    }

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mTag;
    }
};

class ParentClass
{
public:
    unsigned mTag;
    ChildClass* mpChild;
    ParentClass(ChildClass* pChild) : mTag(0), mpChild(pChild)
    {
        mpChild->SetParent(this);
    }

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mTag;
    }
};
#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(ParentClass)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance of ParentClass.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const ParentClass * t, const unsigned int file_version)
{
    ar << t->mpChild;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a ParentClass instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, ParentClass * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of ParentClass.
    ChildClass* p_child;
    ar >> p_child;
    ::new(t)ParentClass(p_child);
}
}
} // namespace ...

class ClassOfSimpleVariables
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mNumber;
        archive & mString;
        archive & mVectorOfDoubles; // include <boost/serialization/vector.hpp> for this
        archive & mVectorOfBools;
    }

    int mNumber;
    std::string mString;
    std::vector<double> mVectorOfDoubles;
    std::vector<bool> mVectorOfBools;

public:

    ClassOfSimpleVariables()
    {
        //Do nothing.  Used when loading into a pointer
    }
    ClassOfSimpleVariables(int initial,
                           std::string string,
                           std::vector<double> doubles,
                           std::vector<bool> bools)
            : mString(string),
            mVectorOfDoubles(doubles),
            mVectorOfBools(bools)
    {
        mNumber = initial;
    }

    int GetNumber() const
    {
        return mNumber;
    }

    std::string GetString()
    {
        return mString;
    }

    std::vector<double>& GetVectorOfDoubles()
    {
        return mVectorOfDoubles;
    }

    std::vector<bool>& GetVectorOfBools()
    {
        return mVectorOfBools;
    }
};

class TestArchiving : public CxxTest::TestSuite
{
public:
    void TestArchiveSimpleVars()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simple_vars.arch";

        // Create an output archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            ClassOfSimpleVariables i(42,"hello",doubles,bools);

            // cast to const.
            output_arch << static_cast<const ClassOfSimpleVariables&>(i);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            std::vector<double> bad_doubles(1);
            bad_doubles[0] = 10.3;

            std::vector<bool> bad_bools(1);
            bad_bools[0] = false;

            ClassOfSimpleVariables j(0,"bye",bad_doubles,bad_bools);

            // read the archive
            input_arch >> j;

            // Check that the values
            TS_ASSERT_EQUALS(j.GetNumber(),42);
            TS_ASSERT_EQUALS(j.GetString(),"hello");
            TS_ASSERT_EQUALS(j.GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(j.GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[1],1.2,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT_EQUALS(j.GetVectorOfBools()[0],true);
            TS_ASSERT_EQUALS(j.GetVectorOfBools()[1],true);
        }
    }

    void TestArchivingLinkedChildAndParent() throw (Exception)
    {
        // This test is an abstraction of archiving a cyclically linked parent-child pair.
        // The parent represents a TissueCell and the child represents an AbstractCellCycleModel

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "linked_classes.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ChildClass* p_child = new ChildClass;
            ParentClass* p_parent = new ParentClass(p_child);

            p_child->mTag = 11;
            p_parent->mTag = 10;

            ParentClass* const p_parent_for_archiving = p_parent;
            //ChildClass* const p_child_for_archiving = p_child;

            //output_arch << p_child_for_archiving;
            output_arch << p_parent_for_archiving;

            delete p_child;
            delete p_parent;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            //ChildClass* p_child;
            ParentClass* p_parent;

            input_arch >> p_parent;

            TS_ASSERT_EQUALS(p_parent->mTag, 10u);
            TS_ASSERT_EQUALS(p_parent->mpChild->mTag, 11u);
            TS_ASSERT_EQUALS(p_parent->mpChild->mpParent, p_parent);

            delete p_parent->mpChild;
            delete p_parent;
        }
    }


    void TestArchivingSetOfSetOfPointers() throw (Exception)
    {
        // This test is an abstraction of archiving a set of sets of pointers and a list of objects.
        //
        // Note that the list.push_back method uses the copy constructor. This is why we iterate
        // through the list to generate the pointers to populate the set.

        std::vector<double> doubles;
        std::vector<bool> bools;

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "pointer_set.arch";

        // Save
        {
            // Create aClassOfSimpleVariablesn output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ClassOfSimpleVariables one(42, "hello", doubles,bools);
            ClassOfSimpleVariables two(256, "goodbye", doubles,bools);
            ClassOfSimpleVariables three(1, "not used in set", doubles,bools);

            std::list<ClassOfSimpleVariables> a_list;
            std::set<ClassOfSimpleVariables*> a_set;
            a_list.push_back(one);
            a_set.insert( &(a_list.back()) );
            a_list.push_back(two);
            a_set.insert( &(a_list.back()) );
            a_list.push_back(three);

            std::set<std::set<ClassOfSimpleVariables*> > wrapper_set;
            wrapper_set.insert(a_set);

            output_arch << static_cast<const std::list<ClassOfSimpleVariables>&>(a_list);
            output_arch << static_cast<const std::set<std::set<ClassOfSimpleVariables*> >&>(wrapper_set);
        }

        //Load
        {
            std::set<std::set<ClassOfSimpleVariables*> > wrapper_set;
            std::list<ClassOfSimpleVariables> a_list;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TS_ASSERT_EQUALS(wrapper_set.size(), 0u);
            input_arch >> a_list;
            input_arch >> wrapper_set;
            TS_ASSERT_EQUALS(wrapper_set.size(), 1u);
            const std::set<ClassOfSimpleVariables*>& a_set = *(wrapper_set.begin());
            TS_ASSERT_EQUALS(a_set.size(), 2u);

            ClassOfSimpleVariables* p_one_in_set = NULL;
            ClassOfSimpleVariables* p_two_in_set = NULL;
            for (std::set<ClassOfSimpleVariables*>::iterator it = a_set.begin();
                 it!=a_set.end();
                 ++it)
            {
                   ClassOfSimpleVariables* p_class = *(it);
                   if (p_class->GetNumber()==42)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 42);
                        TS_ASSERT_EQUALS(p_class->GetString(), "hello");
                        p_one_in_set = p_class;
                   }
                   else
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 256);
                        TS_ASSERT_EQUALS(p_class->GetString(), "goodbye");
                        p_two_in_set = p_class;
                   }
            }

            ClassOfSimpleVariables* p_one_in_list = NULL;
            ClassOfSimpleVariables* p_two_in_list = NULL;
            for (std::list<ClassOfSimpleVariables>::iterator it = a_list.begin();
                 it!=a_list.end();
                 ++it)
            {
                   ClassOfSimpleVariables* p_class = &(*it);
                   if (p_class->GetNumber()==42)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 42);
                        TS_ASSERT_EQUALS(p_class->GetString(), "hello");
                        p_one_in_list=p_class;
                   }
                   else if (p_class->GetNumber()==256)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 256);
                        TS_ASSERT_EQUALS(p_class->GetString(), "goodbye");
                        p_two_in_list=p_class;
                   }
                   else
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 1);
                        TS_ASSERT_EQUALS(p_class->GetString(), "not used in set");
                   }
            }
            TS_ASSERT_DIFFERS(p_one_in_list, (void*)NULL);
            TS_ASSERT_DIFFERS(p_two_in_list, (void*)NULL);
            TS_ASSERT_EQUALS(p_one_in_list, p_one_in_set);
            TS_ASSERT_EQUALS(p_two_in_list, p_two_in_set);
        }
    }
};


#endif /*TESTARCHIVING_HPP_*/

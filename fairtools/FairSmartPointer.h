/********************************************************************************
 *    Copyright (C) 2020 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#ifndef FAIRSMARTPOINTER_H_
#define FAIRSMARTPOINTER_H_

#include <functional>
#include <memory>
#include <type_traits>

namespace fair {

// make_unique implementation, until C++14 is required
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// make_unique implementation (array variant), until C++14 is required
template<typename T>
std::unique_ptr<T> make_unique(std::size_t size)
{
    return std::unique_ptr<T>(new typename std::remove_extent<T>::type[size]());
}

}   // namespace fair

#endif /* FAIRSMARTPOINTER_H_ */

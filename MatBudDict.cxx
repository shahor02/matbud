// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MatBudDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Ray.h"
#include "MatLayerCyl.h"
#include "MatLayerCylSet.h"
#include "Ray.h"
#include "MatLayerCyl.h"
#include "MatLayerCylSet.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_Ray(void *p);
   static void deleteArray_Ray(void *p);
   static void destruct_Ray(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Ray*)
   {
      ::Ray *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Ray >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Ray", ::Ray::Class_Version(), "Ray.h", 36,
                  typeid(::Ray), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Ray::Dictionary, isa_proxy, 4,
                  sizeof(::Ray) );
      instance.SetDelete(&delete_Ray);
      instance.SetDeleteArray(&deleteArray_Ray);
      instance.SetDestructor(&destruct_Ray);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Ray*)
   {
      return GenerateInitInstanceLocal((::Ray*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Ray*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MatLayerCyl(void *p = 0);
   static void *newArray_MatLayerCyl(Long_t size, void *p);
   static void delete_MatLayerCyl(void *p);
   static void deleteArray_MatLayerCyl(void *p);
   static void destruct_MatLayerCyl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MatLayerCyl*)
   {
      ::MatLayerCyl *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MatLayerCyl >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MatLayerCyl", ::MatLayerCyl::Class_Version(), "MatLayerCyl.h", 64,
                  typeid(::MatLayerCyl), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MatLayerCyl::Dictionary, isa_proxy, 4,
                  sizeof(::MatLayerCyl) );
      instance.SetNew(&new_MatLayerCyl);
      instance.SetNewArray(&newArray_MatLayerCyl);
      instance.SetDelete(&delete_MatLayerCyl);
      instance.SetDeleteArray(&deleteArray_MatLayerCyl);
      instance.SetDestructor(&destruct_MatLayerCyl);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MatLayerCyl*)
   {
      return GenerateInitInstanceLocal((::MatLayerCyl*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MatLayerCyl*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MatLayerCylSet(void *p = 0);
   static void *newArray_MatLayerCylSet(Long_t size, void *p);
   static void delete_MatLayerCylSet(void *p);
   static void deleteArray_MatLayerCylSet(void *p);
   static void destruct_MatLayerCylSet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MatLayerCylSet*)
   {
      ::MatLayerCylSet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MatLayerCylSet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MatLayerCylSet", ::MatLayerCylSet::Class_Version(), "MatLayerCylSet.h", 27,
                  typeid(::MatLayerCylSet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MatLayerCylSet::Dictionary, isa_proxy, 4,
                  sizeof(::MatLayerCylSet) );
      instance.SetNew(&new_MatLayerCylSet);
      instance.SetNewArray(&newArray_MatLayerCylSet);
      instance.SetDelete(&delete_MatLayerCylSet);
      instance.SetDeleteArray(&deleteArray_MatLayerCylSet);
      instance.SetDestructor(&destruct_MatLayerCylSet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MatLayerCylSet*)
   {
      return GenerateInitInstanceLocal((::MatLayerCylSet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MatLayerCylSet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Ray::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Ray::Class_Name()
{
   return "Ray";
}

//______________________________________________________________________________
const char *Ray::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Ray*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Ray::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Ray*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Ray::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Ray*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Ray::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Ray*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MatLayerCyl::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MatLayerCyl::Class_Name()
{
   return "MatLayerCyl";
}

//______________________________________________________________________________
const char *MatLayerCyl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCyl*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MatLayerCyl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCyl*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MatLayerCyl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCyl*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MatLayerCyl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCyl*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MatLayerCylSet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MatLayerCylSet::Class_Name()
{
   return "MatLayerCylSet";
}

//______________________________________________________________________________
const char *MatLayerCylSet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCylSet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MatLayerCylSet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCylSet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MatLayerCylSet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCylSet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MatLayerCylSet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MatLayerCylSet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Ray::Streamer(TBuffer &R__b)
{
   // Stream an object of class Ray.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Ray::Class(),this);
   } else {
      R__b.WriteClassBuffer(Ray::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Ray(void *p) {
      delete ((::Ray*)p);
   }
   static void deleteArray_Ray(void *p) {
      delete [] ((::Ray*)p);
   }
   static void destruct_Ray(void *p) {
      typedef ::Ray current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Ray

//______________________________________________________________________________
void MatLayerCyl::Streamer(TBuffer &R__b)
{
   // Stream an object of class MatLayerCyl.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MatLayerCyl::Class(),this);
   } else {
      R__b.WriteClassBuffer(MatLayerCyl::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MatLayerCyl(void *p) {
      return  p ? new(p) ::MatLayerCyl : new ::MatLayerCyl;
   }
   static void *newArray_MatLayerCyl(Long_t nElements, void *p) {
      return p ? new(p) ::MatLayerCyl[nElements] : new ::MatLayerCyl[nElements];
   }
   // Wrapper around operator delete
   static void delete_MatLayerCyl(void *p) {
      delete ((::MatLayerCyl*)p);
   }
   static void deleteArray_MatLayerCyl(void *p) {
      delete [] ((::MatLayerCyl*)p);
   }
   static void destruct_MatLayerCyl(void *p) {
      typedef ::MatLayerCyl current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MatLayerCyl

//______________________________________________________________________________
void MatLayerCylSet::Streamer(TBuffer &R__b)
{
   // Stream an object of class MatLayerCylSet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MatLayerCylSet::Class(),this);
   } else {
      R__b.WriteClassBuffer(MatLayerCylSet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MatLayerCylSet(void *p) {
      return  p ? new(p) ::MatLayerCylSet : new ::MatLayerCylSet;
   }
   static void *newArray_MatLayerCylSet(Long_t nElements, void *p) {
      return p ? new(p) ::MatLayerCylSet[nElements] : new ::MatLayerCylSet[nElements];
   }
   // Wrapper around operator delete
   static void delete_MatLayerCylSet(void *p) {
      delete ((::MatLayerCylSet*)p);
   }
   static void deleteArray_MatLayerCylSet(void *p) {
      delete [] ((::MatLayerCylSet*)p);
   }
   static void destruct_MatLayerCylSet(void *p) {
      typedef ::MatLayerCylSet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MatLayerCylSet

namespace ROOT {
   static TClass *vectorlEshortgR_Dictionary();
   static void vectorlEshortgR_TClassManip(TClass*);
   static void *new_vectorlEshortgR(void *p = 0);
   static void *newArray_vectorlEshortgR(Long_t size, void *p);
   static void delete_vectorlEshortgR(void *p);
   static void deleteArray_vectorlEshortgR(void *p);
   static void destruct_vectorlEshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<short>*)
   {
      vector<short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<short>", -2, "vector", 216,
                  typeid(vector<short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshortgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<short>) );
      instance.SetNew(&new_vectorlEshortgR);
      instance.SetNewArray(&newArray_vectorlEshortgR);
      instance.SetDelete(&delete_vectorlEshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshortgR);
      instance.SetDestructor(&destruct_vectorlEshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<short> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<short>*)0x0)->GetClass();
      vectorlEshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<short> : new vector<short>;
   }
   static void *newArray_vectorlEshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<short>[nElements] : new vector<short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshortgR(void *p) {
      delete ((vector<short>*)p);
   }
   static void deleteArray_vectorlEshortgR(void *p) {
      delete [] ((vector<short>*)p);
   }
   static void destruct_vectorlEshortgR(void *p) {
      typedef vector<short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<short>

namespace ROOT {
   static TClass *vectorlEMatLayerCylgR_Dictionary();
   static void vectorlEMatLayerCylgR_TClassManip(TClass*);
   static void *new_vectorlEMatLayerCylgR(void *p = 0);
   static void *newArray_vectorlEMatLayerCylgR(Long_t size, void *p);
   static void delete_vectorlEMatLayerCylgR(void *p);
   static void deleteArray_vectorlEMatLayerCylgR(void *p);
   static void destruct_vectorlEMatLayerCylgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MatLayerCyl>*)
   {
      vector<MatLayerCyl> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MatLayerCyl>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MatLayerCyl>", -2, "vector", 216,
                  typeid(vector<MatLayerCyl>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMatLayerCylgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<MatLayerCyl>) );
      instance.SetNew(&new_vectorlEMatLayerCylgR);
      instance.SetNewArray(&newArray_vectorlEMatLayerCylgR);
      instance.SetDelete(&delete_vectorlEMatLayerCylgR);
      instance.SetDeleteArray(&deleteArray_vectorlEMatLayerCylgR);
      instance.SetDestructor(&destruct_vectorlEMatLayerCylgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MatLayerCyl> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MatLayerCyl>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMatLayerCylgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MatLayerCyl>*)0x0)->GetClass();
      vectorlEMatLayerCylgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMatLayerCylgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMatLayerCylgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatLayerCyl> : new vector<MatLayerCyl>;
   }
   static void *newArray_vectorlEMatLayerCylgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatLayerCyl>[nElements] : new vector<MatLayerCyl>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMatLayerCylgR(void *p) {
      delete ((vector<MatLayerCyl>*)p);
   }
   static void deleteArray_vectorlEMatLayerCylgR(void *p) {
      delete [] ((vector<MatLayerCyl>*)p);
   }
   static void destruct_vectorlEMatLayerCylgR(void *p) {
      typedef vector<MatLayerCyl> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MatLayerCyl>

namespace ROOT {
   static TClass *vectorlEMatCellgR_Dictionary();
   static void vectorlEMatCellgR_TClassManip(TClass*);
   static void *new_vectorlEMatCellgR(void *p = 0);
   static void *newArray_vectorlEMatCellgR(Long_t size, void *p);
   static void delete_vectorlEMatCellgR(void *p);
   static void deleteArray_vectorlEMatCellgR(void *p);
   static void destruct_vectorlEMatCellgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MatCell>*)
   {
      vector<MatCell> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MatCell>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MatCell>", -2, "vector", 216,
                  typeid(vector<MatCell>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMatCellgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<MatCell>) );
      instance.SetNew(&new_vectorlEMatCellgR);
      instance.SetNewArray(&newArray_vectorlEMatCellgR);
      instance.SetDelete(&delete_vectorlEMatCellgR);
      instance.SetDeleteArray(&deleteArray_vectorlEMatCellgR);
      instance.SetDestructor(&destruct_vectorlEMatCellgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MatCell> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MatCell>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMatCellgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MatCell>*)0x0)->GetClass();
      vectorlEMatCellgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMatCellgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMatCellgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatCell> : new vector<MatCell>;
   }
   static void *newArray_vectorlEMatCellgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatCell>[nElements] : new vector<MatCell>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMatCellgR(void *p) {
      delete ((vector<MatCell>*)p);
   }
   static void deleteArray_vectorlEMatCellgR(void *p) {
      delete [] ((vector<MatCell>*)p);
   }
   static void destruct_vectorlEMatCellgR(void *p) {
      typedef vector<MatCell> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MatCell>

namespace {
  void TriggerDictionaryInitialization_MatBudDict_Impl() {
    static const char* headers[] = {
"Ray.h",
"MatLayerCyl.h",
"MatLayerCylSet.h",
"Ray.h",
"MatLayerCyl.h",
"MatLayerCylSet.h",
0
    };
    static const char* includePaths[] = {
"/home/shahoian/alice/sw/ubuntu1710_x86-64/ROOT/v6-12-04-3/include",
"/home/shahoian/alice/sw/ubuntu1710_x86-64/O2/dev-1/include",
"/home/shahoian/alice/sw/ubuntu1710_x86-64/FairRoot/c210f4bca8-2/include",
"./",
"/home/shahoian/alice/sw/ubuntu1710_x86-64/ROOT/v6-12-04-3/include",
"/home/shahoian/dev/matbud/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MatBudDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$Ray.h")))  Ray;
class __attribute__((annotate("$clingAutoload$MatLayerCyl.h")))  MatLayerCyl;
class __attribute__((annotate("$clingAutoload$MatLayerCylSet.h")))  MatLayerCylSet;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MatBudDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Ray.h"
#include "MatLayerCyl.h"
#include "MatLayerCylSet.h"
#include "Ray.h"
#include "MatLayerCyl.h"
#include "MatLayerCylSet.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"MatLayerCyl", payloadCode, "@",
"MatLayerCylSet", payloadCode, "@",
"Ray", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MatBudDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MatBudDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MatBudDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MatBudDict() {
  TriggerDictionaryInitialization_MatBudDict_Impl();
}

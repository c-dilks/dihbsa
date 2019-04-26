// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DihBsaDict

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
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_DIS(void *p = 0);
   static void *newArray_DIS(Long_t size, void *p);
   static void delete_DIS(void *p);
   static void deleteArray_DIS(void *p);
   static void destruct_DIS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DIS*)
   {
      ::DIS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DIS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DIS", ::DIS::Class_Version(), "DIS.h", 24,
                  typeid(::DIS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DIS::Dictionary, isa_proxy, 4,
                  sizeof(::DIS) );
      instance.SetNew(&new_DIS);
      instance.SetNewArray(&newArray_DIS);
      instance.SetDelete(&delete_DIS);
      instance.SetDeleteArray(&deleteArray_DIS);
      instance.SetDestructor(&destruct_DIS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DIS*)
   {
      return GenerateInitInstanceLocal((::DIS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DIS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_Trajectory(void *p);
   static void deleteArray_Trajectory(void *p);
   static void destruct_Trajectory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Trajectory*)
   {
      ::Trajectory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Trajectory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Trajectory", ::Trajectory::Class_Version(), "Trajectory.h", 25,
                  typeid(::Trajectory), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Trajectory::Dictionary, isa_proxy, 4,
                  sizeof(::Trajectory) );
      instance.SetDelete(&delete_Trajectory);
      instance.SetDeleteArray(&deleteArray_Trajectory);
      instance.SetDestructor(&destruct_Trajectory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Trajectory*)
   {
      return GenerateInitInstanceLocal((::Trajectory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Trajectory*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Dihadron(void *p = 0);
   static void *newArray_Dihadron(Long_t size, void *p);
   static void delete_Dihadron(void *p);
   static void deleteArray_Dihadron(void *p);
   static void destruct_Dihadron(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Dihadron*)
   {
      ::Dihadron *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Dihadron >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Dihadron", ::Dihadron::Class_Version(), "Dihadron.h", 29,
                  typeid(::Dihadron), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Dihadron::Dictionary, isa_proxy, 4,
                  sizeof(::Dihadron) );
      instance.SetNew(&new_Dihadron);
      instance.SetNewArray(&newArray_Dihadron);
      instance.SetDelete(&delete_Dihadron);
      instance.SetDeleteArray(&deleteArray_Dihadron);
      instance.SetDestructor(&destruct_Dihadron);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Dihadron*)
   {
      return GenerateInitInstanceLocal((::Dihadron*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Dihadron*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_EventTree(void *p);
   static void deleteArray_EventTree(void *p);
   static void destruct_EventTree(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventTree*)
   {
      ::EventTree *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventTree >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EventTree", ::EventTree::Class_Version(), "EventTree.h", 30,
                  typeid(::EventTree), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventTree::Dictionary, isa_proxy, 4,
                  sizeof(::EventTree) );
      instance.SetDelete(&delete_EventTree);
      instance.SetDeleteArray(&deleteArray_EventTree);
      instance.SetDestructor(&destruct_EventTree);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventTree*)
   {
      return GenerateInitInstanceLocal((::EventTree*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EventTree*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr DIS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DIS::Class_Name()
{
   return "DIS";
}

//______________________________________________________________________________
const char *DIS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DIS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DIS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DIS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DIS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DIS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DIS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DIS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Trajectory::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Trajectory::Class_Name()
{
   return "Trajectory";
}

//______________________________________________________________________________
const char *Trajectory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Trajectory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Trajectory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Trajectory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Trajectory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Trajectory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Trajectory::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Trajectory*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Dihadron::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Dihadron::Class_Name()
{
   return "Dihadron";
}

//______________________________________________________________________________
const char *Dihadron::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Dihadron*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Dihadron::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Dihadron*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Dihadron::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Dihadron*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Dihadron::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Dihadron*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EventTree::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventTree::Class_Name()
{
   return "EventTree";
}

//______________________________________________________________________________
const char *EventTree::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventTree*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventTree::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventTree*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventTree::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventTree*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventTree::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventTree*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void DIS::Streamer(TBuffer &R__b)
{
   // Stream an object of class DIS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(DIS::Class(),this);
   } else {
      R__b.WriteClassBuffer(DIS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DIS(void *p) {
      return  p ? new(p) ::DIS : new ::DIS;
   }
   static void *newArray_DIS(Long_t nElements, void *p) {
      return p ? new(p) ::DIS[nElements] : new ::DIS[nElements];
   }
   // Wrapper around operator delete
   static void delete_DIS(void *p) {
      delete ((::DIS*)p);
   }
   static void deleteArray_DIS(void *p) {
      delete [] ((::DIS*)p);
   }
   static void destruct_DIS(void *p) {
      typedef ::DIS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DIS

//______________________________________________________________________________
void Trajectory::Streamer(TBuffer &R__b)
{
   // Stream an object of class Trajectory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Trajectory::Class(),this);
   } else {
      R__b.WriteClassBuffer(Trajectory::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Trajectory(void *p) {
      delete ((::Trajectory*)p);
   }
   static void deleteArray_Trajectory(void *p) {
      delete [] ((::Trajectory*)p);
   }
   static void destruct_Trajectory(void *p) {
      typedef ::Trajectory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Trajectory

//______________________________________________________________________________
void Dihadron::Streamer(TBuffer &R__b)
{
   // Stream an object of class Dihadron.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Dihadron::Class(),this);
   } else {
      R__b.WriteClassBuffer(Dihadron::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Dihadron(void *p) {
      return  p ? new(p) ::Dihadron : new ::Dihadron;
   }
   static void *newArray_Dihadron(Long_t nElements, void *p) {
      return p ? new(p) ::Dihadron[nElements] : new ::Dihadron[nElements];
   }
   // Wrapper around operator delete
   static void delete_Dihadron(void *p) {
      delete ((::Dihadron*)p);
   }
   static void deleteArray_Dihadron(void *p) {
      delete [] ((::Dihadron*)p);
   }
   static void destruct_Dihadron(void *p) {
      typedef ::Dihadron current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Dihadron

//______________________________________________________________________________
void EventTree::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventTree.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EventTree::Class(),this);
   } else {
      R__b.WriteClassBuffer(EventTree::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EventTree(void *p) {
      delete ((::EventTree*)p);
   }
   static void deleteArray_EventTree(void *p) {
      delete [] ((::EventTree*)p);
   }
   static void destruct_EventTree(void *p) {
      typedef ::EventTree current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EventTree

namespace {
  void TriggerDictionaryInitialization_DihBsaDict_Impl() {
    static const char* headers[] = {
"Constants.h",
"DIS.h",
"Trajectory.h",
"Dihadron.h",
"EventTree.h",
0
    };
    static const char* includePaths[] = {
"/home/dilks/builds/root/include",
"/home/dilks/j/dihbsa/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DihBsaDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$DIS.h")))  DIS;
class __attribute__((annotate("$clingAutoload$Trajectory.h")))  Trajectory;
class __attribute__((annotate("$clingAutoload$Dihadron.h")))  Dihadron;
class __attribute__((annotate("$clingAutoload$EventTree.h")))  EventTree;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DihBsaDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DIS", payloadCode, "@",
"Dihadron", payloadCode, "@",
"EventTree", payloadCode, "@",
"Trajectory", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DihBsaDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DihBsaDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DihBsaDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DihBsaDict() {
  TriggerDictionaryInitialization_DihBsaDict_Impl();
}

// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME ATOFLibDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "ATOFLib.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TOFComponent(void *p = nullptr);
   static void *newArray_TOFComponent(Long_t size, void *p);
   static void delete_TOFComponent(void *p);
   static void deleteArray_TOFComponent(void *p);
   static void destruct_TOFComponent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TOFComponent*)
   {
      ::TOFComponent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TOFComponent >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TOFComponent", ::TOFComponent::Class_Version(), "ATOFLib.hh", 33,
                  typeid(::TOFComponent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TOFComponent::Dictionary, isa_proxy, 4,
                  sizeof(::TOFComponent) );
      instance.SetNew(&new_TOFComponent);
      instance.SetNewArray(&newArray_TOFComponent);
      instance.SetDelete(&delete_TOFComponent);
      instance.SetDeleteArray(&deleteArray_TOFComponent);
      instance.SetDestructor(&destruct_TOFComponent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TOFComponent*)
   {
      return GenerateInitInstanceLocal(static_cast<::TOFComponent*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TOFComponent*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TOFWindow(void *p = nullptr);
   static void *newArray_TOFWindow(Long_t size, void *p);
   static void delete_TOFWindow(void *p);
   static void deleteArray_TOFWindow(void *p);
   static void destruct_TOFWindow(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TOFWindow*)
   {
      ::TOFWindow *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TOFWindow >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TOFWindow", ::TOFWindow::Class_Version(), "ATOFLib.hh", 59,
                  typeid(::TOFWindow), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TOFWindow::Dictionary, isa_proxy, 4,
                  sizeof(::TOFWindow) );
      instance.SetNew(&new_TOFWindow);
      instance.SetNewArray(&newArray_TOFWindow);
      instance.SetDelete(&delete_TOFWindow);
      instance.SetDeleteArray(&deleteArray_TOFWindow);
      instance.SetDestructor(&destruct_TOFWindow);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TOFWindow*)
   {
      return GenerateInitInstanceLocal(static_cast<::TOFWindow*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TOFWindow*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ATOFProcess(void *p = nullptr);
   static void *newArray_ATOFProcess(Long_t size, void *p);
   static void delete_ATOFProcess(void *p);
   static void deleteArray_ATOFProcess(void *p);
   static void destruct_ATOFProcess(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ATOFProcess*)
   {
      ::ATOFProcess *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ATOFProcess >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ATOFProcess", ::ATOFProcess::Class_Version(), "ATOFLib.hh", 87,
                  typeid(::ATOFProcess), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ATOFProcess::Dictionary, isa_proxy, 4,
                  sizeof(::ATOFProcess) );
      instance.SetNew(&new_ATOFProcess);
      instance.SetNewArray(&newArray_ATOFProcess);
      instance.SetDelete(&delete_ATOFProcess);
      instance.SetDeleteArray(&deleteArray_ATOFProcess);
      instance.SetDestructor(&destruct_ATOFProcess);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ATOFProcess*)
   {
      return GenerateInitInstanceLocal(static_cast<::ATOFProcess*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ATOFProcess*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_GUIclass(void *p);
   static void deleteArray_GUIclass(void *p);
   static void destruct_GUIclass(void *p);
   static void streamer_GUIclass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GUIclass*)
   {
      ::GUIclass *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GUIclass >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("GUIclass", ::GUIclass::Class_Version(), "ATOFLib.hh", 132,
                  typeid(::GUIclass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GUIclass::Dictionary, isa_proxy, 16,
                  sizeof(::GUIclass) );
      instance.SetDelete(&delete_GUIclass);
      instance.SetDeleteArray(&deleteArray_GUIclass);
      instance.SetDestructor(&destruct_GUIclass);
      instance.SetStreamerFunc(&streamer_GUIclass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GUIclass*)
   {
      return GenerateInitInstanceLocal(static_cast<::GUIclass*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::GUIclass*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TOFComponent::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TOFComponent::Class_Name()
{
   return "TOFComponent";
}

//______________________________________________________________________________
const char *TOFComponent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOFComponent*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TOFComponent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOFComponent*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TOFComponent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOFComponent*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TOFComponent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOFComponent*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TOFWindow::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TOFWindow::Class_Name()
{
   return "TOFWindow";
}

//______________________________________________________________________________
const char *TOFWindow::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOFWindow*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TOFWindow::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOFWindow*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TOFWindow::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOFWindow*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TOFWindow::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOFWindow*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ATOFProcess::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ATOFProcess::Class_Name()
{
   return "ATOFProcess";
}

//______________________________________________________________________________
const char *ATOFProcess::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ATOFProcess*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ATOFProcess::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ATOFProcess*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ATOFProcess::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ATOFProcess*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ATOFProcess::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ATOFProcess*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GUIclass::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *GUIclass::Class_Name()
{
   return "GUIclass";
}

//______________________________________________________________________________
const char *GUIclass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GUIclass*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int GUIclass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GUIclass*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GUIclass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GUIclass*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GUIclass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GUIclass*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TOFComponent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TOFComponent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TOFComponent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TOFComponent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TOFComponent(void *p) {
      return  p ? new(p) ::TOFComponent : new ::TOFComponent;
   }
   static void *newArray_TOFComponent(Long_t nElements, void *p) {
      return p ? new(p) ::TOFComponent[nElements] : new ::TOFComponent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TOFComponent(void *p) {
      delete (static_cast<::TOFComponent*>(p));
   }
   static void deleteArray_TOFComponent(void *p) {
      delete [] (static_cast<::TOFComponent*>(p));
   }
   static void destruct_TOFComponent(void *p) {
      typedef ::TOFComponent current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TOFComponent

//______________________________________________________________________________
void TOFWindow::Streamer(TBuffer &R__b)
{
   // Stream an object of class TOFWindow.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TOFWindow::Class(),this);
   } else {
      R__b.WriteClassBuffer(TOFWindow::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TOFWindow(void *p) {
      return  p ? new(p) ::TOFWindow : new ::TOFWindow;
   }
   static void *newArray_TOFWindow(Long_t nElements, void *p) {
      return p ? new(p) ::TOFWindow[nElements] : new ::TOFWindow[nElements];
   }
   // Wrapper around operator delete
   static void delete_TOFWindow(void *p) {
      delete (static_cast<::TOFWindow*>(p));
   }
   static void deleteArray_TOFWindow(void *p) {
      delete [] (static_cast<::TOFWindow*>(p));
   }
   static void destruct_TOFWindow(void *p) {
      typedef ::TOFWindow current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TOFWindow

//______________________________________________________________________________
void ATOFProcess::Streamer(TBuffer &R__b)
{
   // Stream an object of class ATOFProcess.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ATOFProcess::Class(),this);
   } else {
      R__b.WriteClassBuffer(ATOFProcess::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ATOFProcess(void *p) {
      return  p ? new(p) ::ATOFProcess : new ::ATOFProcess;
   }
   static void *newArray_ATOFProcess(Long_t nElements, void *p) {
      return p ? new(p) ::ATOFProcess[nElements] : new ::ATOFProcess[nElements];
   }
   // Wrapper around operator delete
   static void delete_ATOFProcess(void *p) {
      delete (static_cast<::ATOFProcess*>(p));
   }
   static void deleteArray_ATOFProcess(void *p) {
      delete [] (static_cast<::ATOFProcess*>(p));
   }
   static void destruct_ATOFProcess(void *p) {
      typedef ::ATOFProcess current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ATOFProcess

//______________________________________________________________________________
void GUIclass::Streamer(TBuffer &R__b)
{
   // Stream an object of class GUIclass.

   TGMainFrame::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_GUIclass(void *p) {
      delete (static_cast<::GUIclass*>(p));
   }
   static void deleteArray_GUIclass(void *p) {
      delete [] (static_cast<::GUIclass*>(p));
   }
   static void destruct_GUIclass(void *p) {
      typedef ::GUIclass current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_GUIclass(TBuffer &buf, void *obj) {
      ((::GUIclass*)obj)->::GUIclass::Streamer(buf);
   }
} // end of namespace ROOT for class ::GUIclass

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETOFWindowgR_Dictionary();
   static void vectorlETOFWindowgR_TClassManip(TClass*);
   static void *new_vectorlETOFWindowgR(void *p = nullptr);
   static void *newArray_vectorlETOFWindowgR(Long_t size, void *p);
   static void delete_vectorlETOFWindowgR(void *p);
   static void deleteArray_vectorlETOFWindowgR(void *p);
   static void destruct_vectorlETOFWindowgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TOFWindow>*)
   {
      vector<TOFWindow> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TOFWindow>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TOFWindow>", -2, "vector", 389,
                  typeid(vector<TOFWindow>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETOFWindowgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TOFWindow>) );
      instance.SetNew(&new_vectorlETOFWindowgR);
      instance.SetNewArray(&newArray_vectorlETOFWindowgR);
      instance.SetDelete(&delete_vectorlETOFWindowgR);
      instance.SetDeleteArray(&deleteArray_vectorlETOFWindowgR);
      instance.SetDestructor(&destruct_vectorlETOFWindowgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TOFWindow> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TOFWindow>","std::vector<TOFWindow, std::allocator<TOFWindow> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TOFWindow>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETOFWindowgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TOFWindow>*>(nullptr))->GetClass();
      vectorlETOFWindowgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETOFWindowgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETOFWindowgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TOFWindow> : new vector<TOFWindow>;
   }
   static void *newArray_vectorlETOFWindowgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TOFWindow>[nElements] : new vector<TOFWindow>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETOFWindowgR(void *p) {
      delete (static_cast<vector<TOFWindow>*>(p));
   }
   static void deleteArray_vectorlETOFWindowgR(void *p) {
      delete [] (static_cast<vector<TOFWindow>*>(p));
   }
   static void destruct_vectorlETOFWindowgR(void *p) {
      typedef vector<TOFWindow> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TOFWindow>

namespace ROOT {
   static TClass *vectorlETOFComponentgR_Dictionary();
   static void vectorlETOFComponentgR_TClassManip(TClass*);
   static void *new_vectorlETOFComponentgR(void *p = nullptr);
   static void *newArray_vectorlETOFComponentgR(Long_t size, void *p);
   static void delete_vectorlETOFComponentgR(void *p);
   static void deleteArray_vectorlETOFComponentgR(void *p);
   static void destruct_vectorlETOFComponentgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TOFComponent>*)
   {
      vector<TOFComponent> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TOFComponent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TOFComponent>", -2, "vector", 389,
                  typeid(vector<TOFComponent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETOFComponentgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TOFComponent>) );
      instance.SetNew(&new_vectorlETOFComponentgR);
      instance.SetNewArray(&newArray_vectorlETOFComponentgR);
      instance.SetDelete(&delete_vectorlETOFComponentgR);
      instance.SetDeleteArray(&deleteArray_vectorlETOFComponentgR);
      instance.SetDestructor(&destruct_vectorlETOFComponentgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TOFComponent> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TOFComponent>","std::vector<TOFComponent, std::allocator<TOFComponent> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TOFComponent>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETOFComponentgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TOFComponent>*>(nullptr))->GetClass();
      vectorlETOFComponentgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETOFComponentgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETOFComponentgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TOFComponent> : new vector<TOFComponent>;
   }
   static void *newArray_vectorlETOFComponentgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TOFComponent>[nElements] : new vector<TOFComponent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETOFComponentgR(void *p) {
      delete (static_cast<vector<TOFComponent>*>(p));
   }
   static void deleteArray_vectorlETOFComponentgR(void *p) {
      delete [] (static_cast<vector<TOFComponent>*>(p));
   }
   static void destruct_vectorlETOFComponentgR(void *p) {
      typedef vector<TOFComponent> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TOFComponent>

namespace {
  void TriggerDictionaryInitialization_ATOFLibDict_Impl() {
    static const char* headers[] = {
"ATOFLib.hh",
nullptr
    };
    static const char* includePaths[] = {
"/home/terawatt/Programs/root/root-install/include/",
"/home/terawatt/Programs/atoflib/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "ATOFLibDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ATOFLib.hh")))  TOFComponent;
class __attribute__((annotate("$clingAutoload$ATOFLib.hh")))  TOFWindow;
class __attribute__((annotate("$clingAutoload$ATOFLib.hh")))  ATOFProcess;
class __attribute__((annotate("$clingAutoload$ATOFLib.hh")))  GUIclass;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "ATOFLibDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "ATOFLib.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ATOFProcess", payloadCode, "@",
"GUIclass", payloadCode, "@",
"TOFComponent", payloadCode, "@",
"TOFWindow", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ATOFLibDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ATOFLibDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ATOFLibDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ATOFLibDict() {
  TriggerDictionaryInitialization_ATOFLibDict_Impl();
}

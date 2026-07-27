//
// This file is derived from
// sqdbcpp (http://code.google.com/p/sqdbcpp/)
// project. This file is licensed under
// New BSD License (The BSD 2-Clause License).
//
// The code is a little bit modified from the
// original version by Masahiro Kasahara.
//

#ifndef SQDB_SQDB_H
#define SQDB_SQDB_H

#include <string>

#include "sqlite3.h"

#ifdef _WIN32
#  include <tchar.h>
#  define SQDB_MAKE_TEXT(x) _TEXT(x)
#  define SQDB_STRLEN _tcslen
#  define SQDB_STRDUP _tcsdup
#else
#  define SQDB_MAKE_TEXT(x) (x) 
#  define SQDB_STRLEN strlen
#  define SQDB_STRDUP strdup
#endif

#if !defined(SQDB_UTF16) && !defined(SQDB_UTF8)
#  ifdef _WIN32
#    if defined(UNICODE) || defined(_UNICODE)
#      define SQDB_UTF16
#    else
#      define SQDB_UTF8
#    endif
#  else
#    define SQDB_UTF8
#  endif
#endif

#ifdef SQDB_UTF8
#  define SQDB_CHAR char 
#  define SQDB_STD_STRING std::string
#endif

#ifdef SQDB_UTF16
#  define SQDB_CHAR TCHAR 
#  define SQDB_STD_STRING std::wstring
#endif

namespace sqdb
{

// Some conversions below are dangerous unless the caller knows what it is
// doing; make them explicit when the compiler is new enough to support it.
#if defined(SQDB_NO_EXPLICIT_CONVERSION)
#  define SQDB_EXPLICIT_CONVERSION
#elif (defined(__cplusplus) && __cplusplus >= 201103L) || (defined(_MSC_VER) && _MSC_VER >= 1800)
#  define SQDB_EXPLICIT_CONVERSION explicit
#else
#  define SQDB_EXPLICIT_CONVERSION
#endif

class Exception
{
public:
  Exception(sqlite3* db);

  Exception(sqlite3* db, int errorCode);

  Exception(const SQDB_CHAR* errorMsg);

  // Exception owns m_errorMsg, so it needs a deep copy (rule of three).
  // Copying used to hand the same pointer to two objects, both of which
  // free()d it (double free) -- e.g. 'catch (sqdb::Exception e)' by value or
  // a 'throw e;' re-throw.
  Exception(const Exception& x);
  Exception& operator=(const Exception& x);

  ~Exception();

  int GetErrorCode() const;

  // Never returns NULL; returns an empty string if no message is available.
  const SQDB_CHAR* GetErrorMsg() const;
private:
  int m_errorCode;
  SQDB_CHAR* m_errorMsg;
};

// NOTE: wrapped in do { } while (0) so that
//   if (cond) CHECK(db, rc); else foo();
// keeps the 'else' attached to the caller's 'if'.  returnCode is evaluated
// exactly once.
#define CHECK(db, returnCode)                                        \
  do {                                                               \
    const int sqdb_check_return_code_ = (returnCode);                \
    if ( sqdb_check_return_code_ != SQLITE_OK )                      \
      throw ::sqdb::Exception(db, sqdb_check_return_code_);          \
  } while (0)

// Intrusive reference counter shared by the handle wrappers below.
//
// The counted resource itself is owned by the derived class, so the derived
// class must release it when the last reference goes away.  For assignment
// use Reassign(), which drops the reference to the old counter and takes a
// reference to x's counter in one step and tells the caller whether it has
// to release the resource it used to share.  A derived operator= must look
// like:
//
//   Derived& Derived::operator=(const Derived& x)
//   {
//     if ( this != &x ) {
//       Resource* const old = m_resource;          // remember
//       const bool releaseOld = Reassign(x);       // re-bind the counter
//       m_resource = x.m_resource;                 // adopt
//       if ( releaseOld ) Release(old);            // free the old one
//     }
//     return *this;
//   }
class RefCount
{
protected:
  RefCount();

  RefCount(const RefCount& x);

  // Drops the reference to the current counter and takes a reference to
  // x's counter.  Returns true iff this object held the *last* reference to
  // the previous counter, i.e. the caller must now release the resource that
  // counter was guarding.  Safe when *this and x already share a counter.
  bool Reassign(const RefCount& x);

  void IncRef();
  unsigned DecRef();

private:
  // Not implemented on purpose: it cannot tell the derived class that the
  // previously shared resource has to be released.  Use Reassign() instead.
  RefCount& operator=(const RefCount& x);

  unsigned* m_refCount;
};

class Blob : public RefCount
{
public:
  Blob(const void* data, int size);

  Blob(const Blob& x);
  Blob& operator=(const Blob& x);

  int GetSize() const;
  const char* GetData() const;

  ~Blob();

private:
  char* m_data;
  int m_size;
};

class Convertor
{
public:
  Convertor(sqlite3* db, sqlite3_stmt* stmt, int field);

  operator int() const;
  operator long long() const;
  operator double() const;
  operator SQDB_STD_STRING() const;
  // DANGEROUS, see GetText() below.  Explicit (when the compiler supports it)
  // so that it cannot be picked up by accident.
  SQDB_EXPLICIT_CONVERSION operator const SQDB_CHAR*() const;
  operator Blob() const;

  int GetInt() const;
  long long GetLongLong() const;
  double GetDouble() const;
  // Returns an empty string when the column is SQL NULL.
  SQDB_STD_STRING GetString() const;
  // !!! DANGER !!!  The returned pointer belongs to SQLite, NOT to the caller.
  // It is invalidated by the next sqlite3_step()/reset()/finalize() on the
  // statement (i.e. by Statement::Next(), Statement::Bind(), or the
  // destruction of the Statement) and even by converting *another* column of
  // the same row to a different type.  Copy it into a SQDB_STD_STRING (see
  // GetString()) if it has to outlive the current row.  May return NULL when
  // the column is SQL NULL.
  const SQDB_CHAR* GetText() const;
  Blob GetBlob() const;

private:
  sqlite3* m_db;
  sqlite3_stmt* m_stmt;
  int m_field;
};

class Statement : public RefCount
{
public:
  Statement(sqlite3* db, sqlite3_stmt* stmt);

  Statement(const Statement& x);
  Statement& operator=(const Statement& x);

  bool Next();
  Convertor GetField(int field) const;

  template<class T>
  void Bind(int i, const T& value)
  {
    if ( *m_needReset )
      Reset();
    DoBind(i, value);
  }

  void BindBlob(int i, const void* value, int n);
  void BindNull(int i);

  ~Statement();

private:
  void DoBind(int i, int value); 
  void DoBind(int i, long long value); 
  void DoBind(int i, double value);
  void DoBind(int i, const SQDB_STD_STRING& value);
  void DoBind(int i, const SQDB_CHAR* value);

  // Bind blob.
  void DoBind(int i, const void* value, int n);

  // Bind null.
  void DoBind(int i);

  // Reset binders so that new values can be bound.
  void Reset();

  sqlite3* m_db;
  sqlite3_stmt* m_stmt;
  // The "a value has been bound and the statement has been stepped, so it must
  // be reset before the next bind" state belongs to the sqlite3_stmt, not to
  // an individual Statement handle, so every copy has to see the same flag.
  bool* m_needReset;
};

class QueryStr
{
public:
  QueryStr();

  // QueryStr owns m_buf, so it needs a deep copy (rule of three); copying it
  // used to double-free the sqlite-allocated buffer.
  QueryStr(const QueryStr& x);
  QueryStr& operator=(const QueryStr& x);

  const SQDB_CHAR* Format(const SQDB_CHAR* fmt, ...);

  const SQDB_CHAR* Get() const;

  ~QueryStr();

private:
  static SQDB_CHAR* CloneBuf(const SQDB_CHAR* buf);
  static void FreeBuf(SQDB_CHAR* buf);

  SQDB_CHAR* m_buf;
};

class Db : public RefCount
{
public:
  Db(const SQDB_CHAR* fileName);

  void BeginTransaction();
  void CommitTransaction();
  void RollbackTransaction();

  bool TableExists(const SQDB_CHAR* tableName);
  Statement Query(const SQDB_CHAR* queryStr);
  inline void Do(const SQDB_CHAR* queryStr) { Query(queryStr).Next(); }
  long long LastId();
  void MakeItFasterAndDangerous();

  Db(const Db& x);
  Db& operator=(const Db& x);

  ~Db();

private:
  sqlite3* m_db;
};

}

#endif


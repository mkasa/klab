//
// This file is derived from
// sqdbcpp (http://code.google.com/p/sqdbcpp/)
// project. This file is licensed under
// New BSD License (The BSD 2-Clause License).
//
#include <cassert>
#include <climits>
#include <cstdio>
#include <cstring>
#include <memory>

#include <stdlib.h>

#include "sqdb.h"

using namespace sqdb;

namespace {

// SQDB_STRDUP(NULL) is undefined behaviour and a NULL message would make
// GetErrorMsg() hand out a NULL pointer, so normalize both here.
SQDB_CHAR* DupErrorMsg(const SQDB_CHAR* msg)
{
  if ( msg == NULL )
    msg = SQDB_MAKE_TEXT("");
  return SQDB_STRDUP(msg);
}

const SQDB_CHAR* SqliteErrorMsg(sqlite3* db)
{
  return (const SQDB_CHAR*)
#ifdef SQDB_UTF8
    sqlite3_errmsg
#else
    sqlite3_errmsg16
#endif
    (db);
}

} // anonymous namespace

Exception::Exception(sqlite3* db)
: m_errorCode(sqlite3_errcode(db)),
  m_errorMsg(DupErrorMsg(SqliteErrorMsg(db)))
{
}

Exception::Exception(sqlite3* db, int errorCode)
: m_errorCode(errorCode),
  m_errorMsg(DupErrorMsg(SqliteErrorMsg(db)))
{
}

Exception::Exception(const SQDB_CHAR* errorMsg)
: m_errorCode(-1),
  m_errorMsg(DupErrorMsg(errorMsg))
{
}

Exception::Exception(const Exception& x)
: m_errorCode(x.m_errorCode),
  m_errorMsg(DupErrorMsg(x.m_errorMsg))
{
}

Exception& Exception::operator=(const Exception& x)
{
  if ( this != &x )
  {
    SQDB_CHAR* const newMsg = DupErrorMsg(x.m_errorMsg);
    free(m_errorMsg);
    m_errorMsg = newMsg;
    m_errorCode = x.m_errorCode;
  }
  return *this;
}

Exception::~Exception()
{
  free(m_errorMsg);
}

int Exception::GetErrorCode() const
{
  return m_errorCode;
}

const SQDB_CHAR* Exception::GetErrorMsg() const
{
  // m_errorMsg is NULL only if strdup() itself failed (out of memory).
  return m_errorMsg != NULL ? m_errorMsg : SQDB_MAKE_TEXT("");
}

RefCount::RefCount()
: m_refCount(NULL)
{
}

RefCount::RefCount(const RefCount& x)
: m_refCount(x.m_refCount)
{
}

bool RefCount::Reassign(const RefCount& x)
{
  if ( this == &x )
    return false;

  // Drop the reference to the counter we are bound to right now.  Note that
  // x may well share that very counter, in which case its count is >= 2 and
  // the decrement below cannot delete it.
  bool wasLastReference = false;
  if ( m_refCount != NULL && --*m_refCount == 0 )
  {
    delete m_refCount;
    wasLastReference = true;
  }
  m_refCount = NULL;

  assert(x.m_refCount != NULL);
  m_refCount = x.m_refCount;
  IncRef();

  return wasLastReference;
}

void RefCount::IncRef()
{
  if ( !m_refCount ) 
  {
    m_refCount = new unsigned;
    *m_refCount = 0;
  }
  ++*m_refCount;
}

unsigned RefCount::DecRef()
{
  assert(m_refCount);
  unsigned value = --*m_refCount;
  if ( value == 0 ) 
  {
    delete m_refCount;
    m_refCount = NULL;
  }
  return value;
}

Blob::Blob(const void* data, int size)
: m_size(size)
{
  m_data = new char[size];
  std::uninitialized_copy((char*)data, (char*)data + size, m_data); 
  IncRef();
}

Blob::Blob(const Blob& x)
: RefCount(x), m_data(x.m_data), m_size(x.m_size)
{
  IncRef();
}

Blob& Blob::operator=(const Blob& x)
{
  if ( this != &x )
  {
    char* const oldData = m_data;
    const bool releaseOld = Reassign(x);
    m_data = x.m_data;
    m_size = x.m_size;
    if ( releaseOld )
      delete[] oldData;
  }
  return *this;
}

int Blob::GetSize() const
{
  return m_size;
}

const char* Blob::GetData() const
{
  return m_data;
}

Blob::~Blob()
{
  if ( DecRef() == 0 ) 
  {
    delete[] m_data;
  }
}

Convertor::Convertor(sqlite3* db, sqlite3_stmt* stmt, int field)
: m_db(db), m_stmt(stmt), m_field(field)
{
}

Convertor::operator int() const
{
  return GetInt();
}

Convertor::operator long long() const
{
  return GetLongLong();
}

Convertor::operator double() const
{
  return GetDouble();
}

Convertor::operator SQDB_STD_STRING() const
{
  return GetString();
}

Convertor::operator const SQDB_CHAR*() const
{
  return GetText();
}

Convertor::operator Blob() const
{
  return GetBlob();
}

int Convertor::GetInt() const
{
  assert(m_stmt);
  return sqlite3_column_int(m_stmt, m_field);
}

long long Convertor::GetLongLong() const
{
  assert(m_stmt);
  return sqlite3_column_int64(m_stmt, m_field);
}

double Convertor::GetDouble() const
{
  assert(m_stmt);
  return sqlite3_column_double(m_stmt, m_field);
}

SQDB_STD_STRING Convertor::GetString() const
{
  assert(m_stmt);
  const SQDB_CHAR* result = (const SQDB_CHAR*)
#ifdef SQDB_UTF8
  sqlite3_column_text
#else
  sqlite3_column_text16
#endif
  (m_stmt, m_field);
  // sqlite3_column_text() returns NULL for a SQL NULL column (and on OOM);
  // constructing a std::string from NULL is undefined behaviour.
  if ( result == NULL )
    return SQDB_STD_STRING();
  return SQDB_STD_STRING(result);
}

// See the loud warning at the declaration in sqdb.h: the returned pointer is
// owned by SQLite and is invalidated by the next step/reset/finalize of the
// statement, or by converting another column of the same row.
const SQDB_CHAR* Convertor::GetText() const
{
  assert(m_stmt);
  const SQDB_CHAR* result = (const SQDB_CHAR*)
#ifdef SQDB_UTF8
  sqlite3_column_text
#else
  sqlite3_column_text16
#endif
  (m_stmt, m_field);
  return result;
}

Blob Convertor::GetBlob() const
{
  assert(m_stmt);
  const void* data = sqlite3_column_blob(m_stmt, m_field);
  int size = sqlite3_column_bytes(m_stmt, m_field);
  return Blob(data, size);
}

Statement::Statement(sqlite3* db, sqlite3_stmt* stmt)
: RefCount(), m_db(db), m_stmt(stmt), m_needReset(NULL)
{
  m_needReset = new bool(false);
  IncRef();
}

// The reset state is a property of the shared sqlite3_stmt, so copies must
// share it, otherwise a copy would happily bind to a statement that has been
// stepped but not reset (sqlite3_bind_*() then returns SQLITE_MISUSE).
Statement::Statement(const Statement& x)
: RefCount(x), m_db(x.m_db), m_stmt(x.m_stmt), m_needReset(x.m_needReset)
{
  IncRef();
}

Statement& Statement::operator=(const Statement& x)
{
  if ( this != &x )
  {
    sqlite3_stmt* const oldStmt = m_stmt;
    bool* const oldNeedReset = m_needReset;
    const bool releaseOld = Reassign(x);
    m_db = x.m_db;
    m_stmt = x.m_stmt;
    m_needReset = x.m_needReset;
    if ( releaseOld )
    {
      sqlite3_finalize(oldStmt);
      delete oldNeedReset;
    }
  }
  return *this;
}

bool Statement::Next()
{
  assert(m_stmt);
  int ret = sqlite3_step(m_stmt);
  *m_needReset = true;
  if ( ret == SQLITE_DONE )
  {
    return false;
  }
  else if ( ret == SQLITE_ROW )
  {
    return true;
  }
  else
  {
    CHECK(m_db, ret);
    // Not reachable: sqlite3_step() never returns SQLITE_OK, so the CHECK
    // above always throws.  Falling off the end of a non-void function would
    // be undefined behaviour, so be explicit about it.
    throw Exception(SQDB_MAKE_TEXT("Unexpected return value from sqlite3_step()"));
  }
}

Convertor Statement::GetField(int field) const
{
  return Convertor(m_db, m_stmt, field);
}

Statement::~Statement()
{
  if ( DecRef() == 0 )
  {
    sqlite3_finalize(m_stmt);
    delete m_needReset;
  }
}

void Statement::BindBlob(int i, const void* value, int n)
{
  if ( *m_needReset )
    Reset();
  DoBind(i, value, n);
}

void Statement::BindNull(int i)
{
  if ( *m_needReset )
    Reset();
  DoBind(i);
}

void Statement::DoBind(int i, int value)
{
  const int ret = sqlite3_bind_int(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, long long value)
{
  const int ret = sqlite3_bind_int64(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, double value)
{
  const int ret = sqlite3_bind_double(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const SQDB_STD_STRING& value)
{
  // sqlite3_bind_text() takes the byte count as an int.  Truncating size_t to
  // int would silently cut long values short, and a negative value would make
  // SQLite treat the argument as NUL-terminated (which changes the semantics
  // for strings containing embedded NULs).  Refuse instead.
  const size_t charSize = sizeof(SQDB_STD_STRING::value_type);
  if ( value.size() > (size_t)INT_MAX / charSize )
    throw Exception(SQDB_MAKE_TEXT("String is too long to be bound (more than INT_MAX bytes)"));
  const int numBytes = (int)(value.size() * charSize);

  const int ret =
#ifdef SQDB_UTF8
  sqlite3_bind_text
#else
  sqlite3_bind_text16
#endif
  (m_stmt, i, value.c_str(), numBytes, SQLITE_TRANSIENT);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const SQDB_CHAR* value)
{
  if ( value == NULL )
  {
    // Binding a NULL pointer as text is meaningless; bind a SQL NULL instead
    // of letting SQDB_STRLEN() dereference it.
    DoBind(i);
    return;
  }

  // See the comment in the SQDB_STD_STRING overload above; SQDB_STRLEN()
  // returns a size_t and the byte count must fit into an int.
  const size_t len = SQDB_STRLEN(value);
  if ( len > (size_t)INT_MAX / sizeof(SQDB_CHAR) )
    throw Exception(SQDB_MAKE_TEXT("String is too long to be bound (more than INT_MAX bytes)"));
  const int numBytes = (int)(len * sizeof(SQDB_CHAR));

  const int ret =
#ifdef SQDB_UTF8
  sqlite3_bind_text
#else
  sqlite3_bind_text16
#endif
  (m_stmt, i, value, numBytes, SQLITE_TRANSIENT);

  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const void* value, int n)
{
  const int ret = sqlite3_bind_blob(m_stmt, i, value, n, SQLITE_TRANSIENT);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i)
{
  const int ret = sqlite3_bind_null(m_stmt, i);
  CHECK(m_db, ret);
}

void Statement::Reset()
{
  assert(m_needReset && *m_needReset);
  assert(m_stmt);

  const int ret = sqlite3_reset(m_stmt);
  // The statement is reset whether or not an error is reported, so clear the
  // flag before (possibly) throwing.
  *m_needReset = false;
  // sqlite3_reset() reports the error code of the previous failed step, which
  // is the only place where that error shows up for some statements; do not
  // throw it away.
  CHECK(m_db, ret);
}

QueryStr::QueryStr()
: m_buf(NULL)
{
}

SQDB_CHAR* QueryStr::CloneBuf(const SQDB_CHAR* buf)
{
  if ( buf == NULL )
    return NULL;
#ifdef SQDB_UTF8
  // m_buf is owned by sqlite3, so the copy has to be allocated by sqlite3 too.
  return sqlite3_mprintf(SQDB_MAKE_TEXT("%s"), buf);
#else
  return SQDB_STRDUP(buf);
#endif
}

void QueryStr::FreeBuf(SQDB_CHAR* buf)
{
#ifdef SQDB_UTF8
  sqlite3_free(buf);
#else
  free(buf);
#endif
}

QueryStr::QueryStr(const QueryStr& x)
: m_buf(CloneBuf(x.m_buf))
{
}

QueryStr& QueryStr::operator=(const QueryStr& x)
{
  if ( this != &x )
  {
    SQDB_CHAR* const newBuf = CloneBuf(x.m_buf);
    FreeBuf(m_buf);
    m_buf = newBuf;
  }
  return *this;
}

const SQDB_CHAR* QueryStr::Format(const SQDB_CHAR* fmt, ...)
{
  va_list va;
  va_start(va, fmt);
#ifdef SQDB_UTF8
  sqlite3_free(m_buf);
  m_buf = sqlite3_vmprintf(fmt, va);
#else
  free(m_buf);
  int len = _vscwprintf(fmt, va) + 1;
  m_buf = (SQDB_CHAR*)malloc(len * sizeof(SQDB_CHAR));
  vswprintf(m_buf, len, fmt, va);
#endif

  va_end(va);

  return m_buf;
}

const SQDB_CHAR* QueryStr::Get() const
{
  return m_buf;
}

QueryStr::~QueryStr()
{
  FreeBuf(m_buf);
}

Db::Db(const SQDB_CHAR* fileName)
: RefCount(), m_db(NULL)
{
#ifdef SQDB_UTF8
  const int ret = sqlite3_open(fileName, &m_db);
#else
  const int ret = sqlite3_open16(fileName, &m_db);
#endif
  if ( ret != SQLITE_OK )
  {
    // sqlite3_open() hands back an allocated handle even when it fails, and
    // this object's destructor will never run because the constructor throws,
    // so the handle (and its file descriptor) has to be closed right here.
    // Build the exception first: it needs the error message from the handle.
    Exception e(m_db, ret);
    sqlite3_close_v2(m_db);
    m_db = NULL;
    throw e;
  }

  IncRef();
}

void Db::BeginTransaction()
{
  Query(SQDB_MAKE_TEXT("BEGIN;")).Next();
}

void Db::CommitTransaction()
{
  Query(SQDB_MAKE_TEXT("COMMIT;")).Next();
}

void Db::RollbackTransaction()
{
  Query(SQDB_MAKE_TEXT("ROLLBACK;")).Next();
}

bool Db::TableExists(const SQDB_CHAR* tableName)
{
  // Bind the name instead of pasting it into the SQL text: sqlite3_vmprintf's
  // "%s" does no quoting at all, so a table name containing an apostrophe
  // produced a syntax error and an attacker-influenced name could inject SQL.
  Statement s =
    Query(SQDB_MAKE_TEXT("select count(*) from sqlite_master where type='table' and name=?;"));
  s.Bind(1, tableName);
  s.Next();
  const int count = s.GetField(0);
  return count > 0;
}

Statement Db::Query(const SQDB_CHAR* queryStr)
{
  if ( queryStr == NULL )
    throw Exception(SQDB_MAKE_TEXT("No SQL statement was given"));

  sqlite3_stmt* stmt = NULL;
  // sqlite3_prepare_v2() (rather than the legacy sqlite3_prepare()) is needed
  // so that sqlite3_step() reports the real error code and message (e.g.
  // "UNIQUE constraint failed: ...") instead of a generic SQLITE_ERROR, and so
  // that statements are re-prepared automatically after a schema change.
#ifdef SQDB_UTF8
  const SQDB_CHAR* tail = NULL;
  const int ret = sqlite3_prepare_v2(m_db, queryStr, -1, &stmt, &tail);
#else
  const void* tail = NULL;
  const int ret = sqlite3_prepare16_v2(m_db, queryStr, -1, &stmt, &tail);
#endif
  CHECK(m_db, ret);

  if ( stmt == NULL )
  {
    // An empty (or comment-only) query string prepares successfully into a
    // NULL statement; stepping it would return SQLITE_MISUSE with an
    // unrelated error message.
    throw Exception(SQDB_MAKE_TEXT("The given query string contains no SQL statement"));
  }

  // Anything after the first statement used to be discarded silently, so
  // Do("stmt1; stmt2;") only ever ran stmt1.  Preparing the remainder is the
  // reliable way to tell "nothing but whitespace/comments left" from "another
  // statement".
  if ( tail != NULL )
  {
    sqlite3_stmt* extraStmt = NULL;
#ifdef SQDB_UTF8
    const int tailRet = sqlite3_prepare_v2(m_db, tail, -1, &extraStmt, NULL);
#else
    const int tailRet = sqlite3_prepare16_v2(m_db, tail, -1, &extraStmt, NULL);
#endif
    if ( tailRet != SQLITE_OK || extraStmt != NULL )
    {
      if ( extraStmt != NULL )
        sqlite3_finalize(extraStmt);
      sqlite3_finalize(stmt);
      if ( tailRet != SQLITE_OK )
        throw Exception(m_db, tailRet);
      throw Exception(SQDB_MAKE_TEXT("Only one SQL statement can be executed at a time; "
                                     "the statement(s) following the first one would be ignored"));
    }
  }

  return Statement(m_db, stmt);
}

long long Db::LastId()
{
  long long ret = sqlite3_last_insert_rowid(m_db);
  return ret;
}

void Db::MakeItFasterAndDangerous()
{
  Do("PRAGMA synchronous = OFF");
  Do("PRAGMA journal_mode = MEMORY");
}

Db::Db(const Db& x)
: RefCount(x),
  m_db(x.m_db)
{
  IncRef();
}

Db& Db::operator=(const Db& x)
{
  if ( this != &x )
  {
    sqlite3* const oldDb = m_db;
    const bool releaseOld = Reassign(x);
    m_db = x.m_db;
    if ( releaseOld )
      sqlite3_close_v2(oldDb);
  }
  return *this;
}

Db::~Db()
{
  if ( DecRef() == 0 )
  {
    // sqlite3_close() returns SQLITE_BUSY and closes *nothing* while any
    // statement prepared on this connection is still alive, which silently
    // leaked the whole connection (and kept the database file locked) when a
    // Statement outlived its Db.  sqlite3_close_v2() instead marks the
    // connection as a zombie and closes it once the last statement is
    // finalized, so the connection is always released eventually.
    sqlite3_close_v2(m_db);
  }
}


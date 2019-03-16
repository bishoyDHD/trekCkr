#ifndef __MIDAS_READER_
#define __MIDAS_READER_

#include "TFile.h"
#include "TTree.h"
#include <ostream>

#include <istream>
#include <fstream>

#ifdef __APPLE__
#include <sys/_types/_int64_t.h>

#else
#include <stdint.h>
#endif

class MRTRunInfo;
class MRTBinaryBlob;
class TBranchElement;


#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
#define ROOTRET int
#else
#define ROOTRET void
#endif


class MIDASindex
{
 public:
  int64_t offset;
  unsigned int size;
};
class MIDAStree;
class MIDASfile:public TFile
{
  friend class MIDAStree;
 public:
  MIDASfile();
  MIDASfile(const char* fname, Option_t* option = "", const char* ftitle = "", Int_t compress = 1);
  virtual ~MIDASfile();
  virtual void	Add(TObject* obj, Bool_t replace = kFALSE);
  virtual void  Append(TObject* obj, Bool_t replace = kFALSE);
  virtual Int_t	AppendKey(TKey* key);
  virtual void	AppendPad(Option_t* option = "");
  virtual void	Browse(TBrowser* b);
  virtual void	Build(TFile* motherFile = 0, TDirectory* motherDir = 0);
  virtual Bool_t cd(const char* path = 0);
  virtual void	Clear(Option_t* option = "");
  virtual TObject*	CloneObject(const TObject* obj, Bool_t autoadd = kTRUE);
  virtual void	Close(Option_t* option = "");
  virtual void	Copy(TObject&) const;
  virtual Bool_t	Cp(const char* dst, Bool_t progressbar = kTRUE, UInt_t buffersize = 1000000);
  virtual TKey*	CreateKey(TDirectory* mother, const TObject* obj, const char* name, Int_t bufsize);
  virtual TKey*	CreateKey(TDirectory* mother, const void* obj, const TClass* cl, const char* name, Int_t bufsize);
  virtual void	Delete(const char* namecycle = "");
  virtual void	DeleteAll(Option_t* option = "");
  virtual void	Draw(Option_t* option = "");
  virtual void	DrawMap(const char* keys = "*", Option_t* option = "");
  virtual void	FillBuffer(char*& buffer);
  virtual TKey*	FindKey(const char* keyname) const;
  virtual TKey*	FindKeyAny(const char* keyname) const;
  virtual TObject*	FindObject(const char* name) const;
  virtual TObject*	FindObject(const TObject* obj) const;
  virtual TObject*	FindObjectAny(const char* name) const;
  virtual TObject*	FindObjectAnyFile(const char* name) const;
  virtual void	Flush();
  virtual TObject*	Get(const char* namecycle);
  virtual Int_t	GetBufferSize() const;
  virtual Long64_t	GetBytesRead() const;
  virtual Long64_t	GetBytesReadExtra() const;
  virtual Int_t	GetBytesToPrefetch() const;
  virtual Long64_t	GetBytesWritten() const;
  virtual TDirectory*	GetDirectory(const char* apath, Bool_t printError = false, const char* funcname = "GetDirectory");
  virtual Long64_t	GetEND() const;    
  virtual const TUrl*	GetEndpointUrl() const;
  virtual Int_t	GetErrno() const;
  virtual TFile*	GetFile() const;
  virtual TKey*	GetKey(const char* name, Short_t cycle = 9999) const;
  virtual TList*	GetList() const;
  virtual TList*	GetListOfKeys() const;
  virtual TObject*	GetMother() const;
  virtual TDirectory*	GetMotherDir() const;
  virtual Int_t	GetNbytesFree() const;
  virtual Int_t	GetNbytesInfo() const;
  virtual Int_t	GetNbytesKeys() const;
  virtual TString	GetNewUrl();
  virtual Int_t	GetNfree() const;
  virtual Int_t	GetNkeys() const;
  virtual Int_t	GetNProcessIDs() const;
  virtual void*	GetObjectChecked(const char* namecycle, const char* classname);
  virtual void*	GetObjectChecked(const char* namecycle, const TClass* cl);
  virtual void*	GetObjectUnchecked(const char* namecycle);
  virtual Option_t*	GetOption() const;
  virtual const char*	GetPath() const;
  virtual const char*	GetPathStatic() const;
  virtual Int_t	GetReadCalls() const;
  virtual Long64_t	GetSeekDir() const;
  virtual Long64_t	GetSeekFree() const;
  virtual Long64_t	GetSeekInfo() const;
  virtual Long64_t	GetSeekKeys() const;
  virtual Long64_t	GetSeekParent() const;
  virtual Long64_t	GetSize() const;
  //virtual TList*	GetStreamerInfoList();
  virtual void	IncrementProcessIDs();
  virtual Bool_t	IsArchive() const;
  virtual Bool_t	IsFolder() const;
  virtual Bool_t	IsModified() const;
  virtual Bool_t	IsOpen() const;
  virtual Bool_t	IsWritable() const;
  virtual void	ls(Option_t* option = "") const;
  virtual void	MakeFree(Long64_t first, Long64_t last);
  virtual void	MakeProject(const char* dirname, const char* classes = "*", Option_t* option = "new");
  virtual void	Map();
  virtual Bool_t	Matches(const char* name);
  virtual TDirectory*	mkdir(const char* name, const char* title = "");
  virtual Bool_t	MustFlush() const;
  virtual TFile*	OpenFile(const char* name, Option_t* option = "", const char* ftitle = "", Int_t compress = 1, Int_t netopt = 0);
  virtual void	Paint(Option_t* option = "");
  virtual void	Print(Option_t* option = "") const;
  virtual void	Purge(Short_t nkeep = 1);
  virtual void	pwd() const;
  virtual void	ReadAll(Option_t* option = "");
  virtual Bool_t	ReadBuffer(char* buf, Int_t len);
  virtual Bool_t	ReadBuffer(char* buf, Long64_t pos, Int_t len);
  virtual Bool_t	ReadBufferAsync(Long64_t offs, Int_t len);
  virtual Bool_t	ReadBuffers(char* buf, Long64_t* pos, Int_t* len, Int_t nbuf);
  virtual void	ReadFree();
  virtual Int_t	ReadKeys(Bool_t forceRead = kTRUE);
  virtual TProcessID*	ReadProcessID(UShort_t pidf);
  virtual void	ReadStreamerInfo();
  virtual Int_t	ReadTObject(TObject* obj, const char* keyname);
  virtual Int_t	Recover();
  virtual void	RecursiveRemove(TObject* obj);
  virtual TObject*	Remove(TObject*);
  virtual Int_t	ReOpen(Option_t* mode);
  virtual void	ResetAfterMerge(TFileMergeInfo*);
  virtual void	ResetErrno() const;
  virtual void	rmdir(const char* name);
  virtual void	Save();
  virtual Int_t	SaveObjectAs(const TObject* obj, const char* filename = "", Option_t* option = "") const;
  virtual void	SaveSelf(Bool_t force = kFALSE);
  virtual void	Seek(Long64_t offset, TFile::ERelativeTo pos = kBeg);

  virtual void	SetBufferSize(Int_t bufsize);
  virtual void	SetCacheRead(TFileCacheRead* cache, TObject* tree = 0, TFile::ECacheAction action = kDisconnect);
  virtual void	SetCacheWrite(TFileCacheWrite* cache);
  virtual void	SetCompressionAlgorithm(Int_t algorithm = 0);
  virtual void	SetCompressionLevel(Int_t level = 1);
  virtual void	SetCompressionSettings(Int_t settings = 1);
  virtual void	SetEND(Long64_t last);
  virtual void	SetModified();
  virtual void	SetMother(TObject* mother);
  virtual void	SetName(const char* newname);
  virtual void	SetOffset(Long64_t offset, TFile::ERelativeTo pos = kBeg);
  virtual void	SetOption(Option_t* option = ">");
  virtual void	SetReadCalls(Int_t readcalls = 0);
  virtual void	SetSeekDir(Long64_t v);
  virtual void	SetTRefAction(TObject* ref, TObject* parent);
  virtual void	SetWritable(Bool_t writable = kTRUE);
  virtual void	ShowStreamerInfo();
  virtual Int_t	Sizeof() const;
  virtual void	UseCache(Int_t maxCacheSize = 10, Int_t pageSize = 0);
  virtual Int_t	Write(const char* name = 0, Int_t opt = 0, Int_t bufsiz = 0);
  virtual Int_t	Write(const char* name = 0, Int_t opt = 0, Int_t bufsiz = 0) const;
  virtual Bool_t	WriteBuffer(const char* buf, Int_t len);
  virtual void	WriteDirHeader();
  virtual void	WriteFree();
  virtual void	WriteHeader();
  virtual void	WriteKeys();
  virtual Int_t	WriteObjectAny(const void* obj, const char* classname, const char* name, Option_t* option = "", Int_t bufsize = 0);
  virtual Int_t	WriteObjectAny(const void* obj, const TClass* cl, const char* name, Option_t* option = "", Int_t bufsize = 0);
  virtual UShort_t	WriteProcessID(TProcessID* pid);
  virtual void	WriteStreamerInfo();
  virtual Int_t	WriteTObject(const TObject* obj, const char* name = 0, Option_t* option = "", Int_t bufsize = 0);


 private:
  std::ifstream * file;

  std::vector< std::map< unsigned int, MIDASindex > > offsetmap;
  MIDAStree *thetree;
  MRTRunInfo *theruninfo;
 public:
 ClassDef(MIDASfile,1);
};





class MIDASevent
{
 public:
  unsigned short id;
  unsigned short mask;

  int serial;

  unsigned int timestamp;
  unsigned int size;
  friend std::ostream & operator<<(std::ostream & os, const MIDASevent& ev);
};


class MIDASbankheader
{
 public:
  unsigned int banksize;
  unsigned int flags;
  friend std::ostream & operator<<(std::ostream & os, const MIDASbankheader& bh);
};


class MIDASbank
{
 public:
  char name[4];
  unsigned int type;
  unsigned int size;
  friend std::ostream & operator<<(std::ostream & os, const MIDASbank& b);
  
};



class MIDAStree:public TTree
{
 public:
  MIDAStree (MIDASfile *f);
  virtual ~MIDAStree ();


  virtual ROOTRET AddBranchToCache(const char* bname, Bool_t subbranches = kFALSE);
  virtual ROOTRET AddBranchToCache(TBranch* branch, Bool_t subbranches = kFALSE);
  virtual TFriendElement*	AddFriend(const char* treename, const char* filename = "");
  virtual TFriendElement*	AddFriend(const char* treename, TFile* file);
  virtual TFriendElement*	AddFriend(TTree* tree, const char* alias = "", Bool_t warn = kFALSE);
  virtual void	AddTotBytes(Int_t tot);
  virtual void	AddZipBytes(Int_t zip);
  virtual Long64_t	AutoSave(Option_t* option = "");
  virtual Int_t	Branch(TList* list, Int_t bufsize = 32000, Int_t splitlevel = 99);
  virtual Int_t	Branch(const char* folder, Int_t bufsize = 32000, Int_t splitlevel = 99);
  virtual Int_t	Branch(TCollection* list, Int_t bufsize = 32000, Int_t splitlevel = 99, const char* name = "");
  virtual TBranch*	Branch(const char* name, void* address, const char* leaflist, Int_t bufsize = 32000);
  virtual TBranch*	BranchOld(const char* name, const char* classname, void* addobj, Int_t bufsize = 32000, Int_t splitlevel = 1);
  virtual TBranch*	BranchRef();
  virtual TBranch*	Bronch(const char* name, const char* classname, void* addobj, Int_t bufsize = 32000, Int_t splitlevel = 99);
  virtual void	Browse(TBrowser*);
  virtual Int_t	BuildIndex(const char* majorname, const char* minorname = "0");
  virtual TFile*	ChangeFile(TFile* file);
  virtual TTree*	CloneTree(Long64_t nentries = -1, Option_t* option = "");
  virtual void	CopyAddresses(TTree*, Bool_t undo = kFALSE);
  virtual Long64_t	CopyEntries(TTree* tree, Long64_t nentries = -1, Option_t* option = "");
  virtual TTree*	CopyTree(const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual TBasket*	CreateBasket(TBranch*);
  virtual void	Delete(Option_t* option = "");
  virtual void	DirectoryAutoAdd(TDirectory*);
  virtual void	Draw(Option_t* opt);
  virtual Long64_t	Draw(const char* varexp, const TCut& selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Long64_t	Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual void	DropBaskets();
  virtual ROOTRET	DropBranchFromCache(const char* bname, Bool_t subbranches = kFALSE);
  virtual ROOTRET	DropBranchFromCache(TBranch* branch, Bool_t subbranches = kFALSE);
  virtual void	DropBuffers(Int_t nbytes);
  virtual Int_t	Fill();
  virtual TBranch*	FindBranch(const char* name);
  virtual TLeaf*	FindLeaf(const char* name);
  virtual Int_t	Fit(const char* funcname, const char* varexp, const char* selection = "", Option_t* option = "", Option_t* goption = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Int_t	FlushBaskets() const;
  virtual const char*	GetAlias(const char* aliasName) const;
  virtual Long64_t	GetAutoFlush() const;
  virtual Long64_t	GetAutoSave() const;
  virtual TBranch*	GetBranch(const char* name);
  virtual TBranchRef*	GetBranchRef() const;
  virtual Bool_t	GetBranchStatus(const char* branchname) const;
  virtual Long64_t	GetCacheSize() const;
  virtual Long64_t	GetChainEntryNumber(Long64_t entry) const;
  virtual Long64_t	GetChainOffset() const;
  virtual TTree::TClusterIterator	GetClusterIterator(Long64_t firstentry);
  virtual Long64_t	GetEntries() const;
  virtual Long64_t	GetEntries(const char* selection);
  virtual Long64_t	GetEntriesFast() const;
  virtual Long64_t	GetEntriesFriend() const;
  virtual Int_t	GetEntry(Long64_t entry = 0, Int_t getall = 0);
  virtual TEntryList*	GetEntryList();
  virtual Long64_t	GetEntryNumber(Long64_t entry) const;
  virtual Long64_t	GetEntryNumberWithBestIndex(Long64_t major, Long64_t minor = 0) const;
  virtual Long64_t	GetEntryNumberWithIndex(Long64_t major, Long64_t minor = 0) const;
  virtual Int_t	GetEntryWithIndex(Int_t major, Int_t minor = 0);
  virtual Long64_t	GetEstimate() const;
  virtual Int_t	GetFileNumber() const;
  virtual TTree*	GetFriend(const char*) const;
  virtual const char*	GetFriendAlias(TTree*) const;
  virtual Int_t*	GetIndex();
  virtual Double_t*	GetIndexValues();
  virtual TIterator*	GetIteratorOnAllLeaves(Bool_t dir = kIterForward);
  virtual TLeaf*	GetLeaf(const char* name);
  virtual TLeaf*	GetLeaf(const char* branchname, const char* leafname);
  virtual TList*	GetListOfAliases() const;
  virtual TObjArray*	GetListOfBranches();
  virtual TList*	GetListOfClones();
  virtual TList*	GetListOfFriends() const;
  virtual TObjArray*	GetListOfLeaves();
  virtual Long64_t	GetMaxEntryLoop() const;
  virtual Double_t	GetMaximum(const char* columname);
  virtual Long64_t	GetMaxVirtualSize() const;
  virtual Double_t	GetMinimum(const char* columname);
  virtual Int_t	GetNbranches();
  virtual Int_t	GetPacketSize() const;
  virtual TVirtualPerfStats*	GetPerfStats() const;
  virtual Long64_t	GetReadEntry() const;
  virtual Long64_t	GetReadEvent() const;
  virtual Int_t	GetScanField() const;
  virtual Long64_t	GetSelectedRows();
  virtual Int_t	GetTimerInterval() const;
  virtual Long64_t	GetTotBytes() const;
  virtual TTree*	GetTree() const;
  virtual TVirtualIndex*	GetTreeIndex() const;
  virtual Int_t	GetTreeNumber() const;
  virtual Int_t	GetUpdate() const;
  virtual TList*	GetUserInfo();
  virtual Double_t*	GetV1();
  virtual Double_t*	GetV2();
  virtual Double_t*	GetV3();
  virtual Double_t*	GetV4();
  virtual Double_t*	GetVal(Int_t i);
  virtual Double_t*	GetW();
  virtual Double_t	GetWeight() const;
  virtual Long64_t	GetZipBytes() const;
  virtual void	IncrementTotalBuffers(Int_t nbytes);
  virtual Bool_t	IsFolder() const;
  virtual Int_t	LoadBaskets(Long64_t maxmemory = 2000000000);
  virtual Long64_t	LoadTree(Long64_t entry);
  virtual Long64_t	LoadTreeFriend(Long64_t entry, TTree* T);
  virtual Int_t	MakeClass(const char* classname = 0, Option_t* option = "");
  virtual Int_t	MakeCode(const char* filename = 0);
  virtual Int_t	MakeProxy(const char* classname, const char* macrofilename = 0, const char* cutfilename = 0, const char* option = 0, Int_t maxUnrolling = 3);
  virtual Int_t	MakeSelector(const char* selector = 0);
  virtual Long64_t	Merge(TCollection* list, Option_t* option = "");
  virtual Long64_t	Merge(TCollection* list, TFileMergeInfo* info);
  virtual Bool_t	Notify();
  virtual void	OptimizeBaskets(ULong64_t maxMemory = 10000000, Float_t minComp = 1.1, Option_t* option = "");
  virtual void	Print(Option_t* option = "") const;
  virtual void	PrintCacheStats(Option_t* option = "") const;
  virtual Long64_t	Process(void* selector, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Long64_t	Process(const char* filename, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Long64_t	Project(const char* hname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual TSQLResult*	Query(const char* varexp = "", const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Long64_t	ReadFile(const char* filename, const char* branchDescriptor = "", char delimiter = ' ');
  virtual Long64_t	ReadStream(std::istream& inputStream, const char* branchDescriptor = "", char delimiter = ' ');
    virtual void	RecursiveRemove(TObject* obj);
  virtual void	Refresh();
  virtual void	RemoveFriend(TTree*);
  virtual void	Reset(Option_t* option = "");
  virtual void	ResetAfterMerge(TFileMergeInfo*);
  virtual void	ResetBranchAddress(TBranch*);
  virtual void	ResetBranchAddresses();
  virtual Long64_t	Scan(const char* varexp = "", const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual Bool_t	SetAlias(const char* aliasName, const char* aliasFormula);
  virtual void	SetAutoFlush(Long64_t autof = -30000000);
  virtual void	SetAutoSave(Long64_t autos = -300000000);
  virtual void	SetBasketSize(const char* bname, Int_t buffsize = 16000);
  Int_t	SetBranchAddress(const char* bname, void** add, TBranch** ptr = 0);
  virtual Int_t	SetBranchAddress(const char* bname, void* add, TClass* realClass, EDataType datatype, Bool_t isptr);
  virtual Int_t	SetBranchAddress(const char* bname, void* add, TBranch** ptr, TClass* realClass, EDataType datatype, Bool_t isptr);
  virtual void	SetBranchStatus(const char* bname, Bool_t status = 1, UInt_t* found = 0);
  virtual ROOTRET	SetCacheEntryRange(Long64_t first, Long64_t last);
  virtual void	SetCacheLearnEntries(Int_t n = 10);
  virtual ROOTRET	SetCacheSize(Long64_t cachesize = -1);
  virtual void	SetChainOffset(Long64_t offset = 0);
  virtual void	SetCircular(Long64_t maxEntries);
  virtual void	SetDebug(Int_t level = 1, Long64_t min = 0, Long64_t max = 9999999);
  virtual void	SetDefaultEntryOffsetLen(Int_t newdefault, Bool_t updateExisting = kFALSE);
  virtual void	SetDirectory(TDirectory* dir);
  virtual Long64_t	SetEntries(Long64_t n = -1);
  virtual void	SetEntryList(TEntryList* list, Option_t* opt = "");
  virtual void	SetEstimate(Long64_t nentries = 1000000);
  virtual void	SetEventList(TEventList* list);
  virtual void	SetFileNumber(Int_t number = 0);
  virtual void	SetMakeClass(Int_t make);
  virtual void	SetMaxEntryLoop(Long64_t maxev = 1000000000);
  virtual void	SetMaxVirtualSize(Long64_t size = 0);
  virtual void	SetName(const char* name);
  virtual void	SetNotify(TObject* obj);
  virtual void	SetObject(const char* name, const char* title);
  virtual void	SetParallelUnzip(Bool_t opt = kTRUE, Float_t RelSize = -1);
  virtual void	SetPerfStats(TVirtualPerfStats* perf);
  virtual void	SetScanField(Int_t n = 50);
  virtual void	SetTimerInterval(Int_t msec = 333);
  virtual void	SetTreeIndex(TVirtualIndex* index);
  virtual void	SetUpdate(Int_t freq = 0);
  virtual void	SetWeight(Double_t w = 1, Option_t* option = "");
  virtual void	Show(Long64_t entry = -1, Int_t lenmax = 20);
  virtual void	StartViewer();

  virtual ROOTRET StopCacheLearningPhase();
  virtual Int_t	UnbinnedFit(const char* funcname, const char* varexp, const char* selection = "", Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  virtual void	UseCurrentStyle();
  virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
  virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;
  
 private:
  MIDASfile *thefile;
  std::map<unsigned int,MRTBinaryBlob*> mappedBranches;
  TBranchElement *dummyelement;
  
 public:
  ClassDef(MIDAStree,1);
};



#endif


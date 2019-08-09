# echo "setup gamConvFacAlg gamConvFacAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtgamConvFacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtgamConvFacAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtgamConvFacAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtgamConvFacAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtgamConvFacAlgtempfile}
  unset cmtgamConvFacAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtgamConvFacAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtgamConvFacAlgtempfile}
unset cmtgamConvFacAlgtempfile
return $cmtsetupstatus


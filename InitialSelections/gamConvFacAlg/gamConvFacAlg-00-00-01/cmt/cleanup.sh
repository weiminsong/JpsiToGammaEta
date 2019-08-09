# echo "cleanup gamConvFacAlg gamConvFacAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtgamConvFacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtgamConvFacAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtgamConvFacAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtgamConvFacAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtgamConvFacAlgtempfile}
  unset cmtgamConvFacAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtgamConvFacAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtgamConvFacAlgtempfile}
unset cmtgamConvFacAlgtempfile
return $cmtcleanupstatus


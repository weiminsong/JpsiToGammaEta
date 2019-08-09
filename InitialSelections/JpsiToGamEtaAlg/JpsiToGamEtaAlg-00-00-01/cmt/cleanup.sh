# echo "cleanup JpsiToGamEtaAlg JpsiToGamEtaAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtJpsiToGamEtaAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtJpsiToGamEtaAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=JpsiToGamEtaAlg -version=JpsiToGamEtaAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtJpsiToGamEtaAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=JpsiToGamEtaAlg -version=JpsiToGamEtaAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtJpsiToGamEtaAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtJpsiToGamEtaAlgtempfile}
  unset cmtJpsiToGamEtaAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtJpsiToGamEtaAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtJpsiToGamEtaAlgtempfile}
unset cmtJpsiToGamEtaAlgtempfile
return $cmtcleanupstatus


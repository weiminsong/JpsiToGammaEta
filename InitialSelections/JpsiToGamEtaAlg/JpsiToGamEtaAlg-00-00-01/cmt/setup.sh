# echo "setup JpsiToGamEtaAlg JpsiToGamEtaAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtJpsiToGamEtaAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtJpsiToGamEtaAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToGamEtaAlg -version=JpsiToGamEtaAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtJpsiToGamEtaAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToGamEtaAlg -version=JpsiToGamEtaAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtJpsiToGamEtaAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtJpsiToGamEtaAlgtempfile}
  unset cmtJpsiToGamEtaAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtJpsiToGamEtaAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtJpsiToGamEtaAlgtempfile}
unset cmtJpsiToGamEtaAlgtempfile
return $cmtsetupstatus


# echo "setup JpsiToGamEta4Alg JpsiToGamEta4Alg-00-00-01 in /besfs/groups/higgs/users/yuanxq/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtJpsiToGamEta4Algtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtJpsiToGamEta4Algtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToGamEta4Alg -version=JpsiToGamEta4Alg-00-00-01 -path=/besfs/groups/higgs/users/yuanxq/workArea/704  -no_cleanup $* >${cmtJpsiToGamEta4Algtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=JpsiToGamEta4Alg -version=JpsiToGamEta4Alg-00-00-01 -path=/besfs/groups/higgs/users/yuanxq/workArea/704  -no_cleanup $* >${cmtJpsiToGamEta4Algtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtJpsiToGamEta4Algtempfile}
  unset cmtJpsiToGamEta4Algtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtJpsiToGamEta4Algtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtJpsiToGamEta4Algtempfile}
unset cmtJpsiToGamEta4Algtempfile
return $cmtsetupstatus


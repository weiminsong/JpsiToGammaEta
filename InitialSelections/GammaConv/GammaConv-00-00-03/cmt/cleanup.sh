# echo "cleanup GammaConv GammaConv-00-00-03 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtGammaConvtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtGammaConvtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=GammaConv -version=GammaConv-00-00-03 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtGammaConvtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=GammaConv -version=GammaConv-00-00-03 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtGammaConvtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtGammaConvtempfile}
  unset cmtGammaConvtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtGammaConvtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtGammaConvtempfile}
unset cmtGammaConvtempfile
return $cmtcleanupstatus


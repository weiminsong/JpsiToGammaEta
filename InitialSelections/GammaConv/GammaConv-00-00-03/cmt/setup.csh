# echo "setup GammaConv GammaConv-00-00-03 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtGammaConvtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtGammaConvtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=GammaConv -version=GammaConv-00-00-03 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtGammaConvtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=GammaConv -version=GammaConv-00-00-03 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtGammaConvtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtGammaConvtempfile}
  unset cmtGammaConvtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtGammaConvtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtGammaConvtempfile}
unset cmtGammaConvtempfile
exit $cmtsetupstatus


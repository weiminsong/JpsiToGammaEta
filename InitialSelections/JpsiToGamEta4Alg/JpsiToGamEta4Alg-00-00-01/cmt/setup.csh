# echo "setup JpsiToGamEta4Alg JpsiToGamEta4Alg-00-00-01 in /besfs/groups/higgs/users/yuanxq/workArea/704"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtJpsiToGamEta4Algtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtJpsiToGamEta4Algtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=JpsiToGamEta4Alg -version=JpsiToGamEta4Alg-00-00-01 -path=/besfs/groups/higgs/users/yuanxq/workArea/704  -no_cleanup $* >${cmtJpsiToGamEta4Algtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=JpsiToGamEta4Alg -version=JpsiToGamEta4Alg-00-00-01 -path=/besfs/groups/higgs/users/yuanxq/workArea/704  -no_cleanup $* >${cmtJpsiToGamEta4Algtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtJpsiToGamEta4Algtempfile}
  unset cmtJpsiToGamEta4Algtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtJpsiToGamEta4Algtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtJpsiToGamEta4Algtempfile}
unset cmtJpsiToGamEta4Algtempfile
exit $cmtsetupstatus


properties([
  disableConcurrentBuilds(),
  pipelineTriggers([
    upstream(
      threshold: 'SUCCESS',
      upstreamProjects: '/TRIQS/triqs/' + env.BRANCH_NAME.replaceAll('/', '%2F')
    )
  ])
])

def platforms = [:]

def dockerPlatforms = ["ubuntu-clang", "ubuntu-gcc", "centos-gcc"]
for (int i = 0; i < dockerPlatforms.size(); i++) {
  def platform = dockerPlatforms[i]
  platforms[platform] = { ->
    stage(platform) {
      timeout(time: 1, unit: 'HOURS') {
	node('docker') {
	  checkout scm
	  /* construct a Dockerfile for this base */
	  sh '''
	    ( echo "FROM flatironinstitute/triqs:$BRANCH_NAME-$STAGE_NAME" ; sed '0,/^FROM /d' Dockerfile ) > Dockerfile.jenkins
	    mv -f Dockerfile.jenkins Dockerfile
	  '''
	  /* build and tag */
	  def img = docker.build("flatironinstitute/dft_tools:${env.BRANCH_NAME}-${env.STAGE_NAME}")
	}
      }
    }
  }
}

parallel platforms

def triqsBranch = env.CHANGE_TARGET ?: env.BRANCH_NAME
def triqsProject = '/TRIQS/triqs/' + triqsBranch.replaceAll('/', '%2F')

properties([
  disableConcurrentBuilds(),
  buildDiscarder(logRotator(numToKeepStr: '10', daysToKeepStr: '30')),
  pipelineTriggers([
    upstream(
      threshold: 'SUCCESS',
      upstreamProjects: triqsProject
    )
  ])
])

def platforms = [:]

def dockerPlatforms = ["ubuntu-clang", "ubuntu-gcc", "centos-gcc"]
for (int i = 0; i < dockerPlatforms.size(); i++) {
  def platform = dockerPlatforms[i]
  platforms[platform] = { -> stage(platform) {
    timeout(time: 1, unit: 'HOURS') {
      node('docker') {
	docker.image("flatironinstitute/triqs:${triqsBranch}-${env.STAGE_NAME}").inside {
	  def srcDir = pwd()

	  checkout scm

	  dir("$srcDir/build") {
	    sh "cmake $srcDir -DTRIQS_ROOT=\$INSTALL -DCMAKE_INSTALL_PREFIX=$srcDir/install"
	    sh "make -j2"
	    try {
	      sh "make test && false"
	    } catch (exc) {
	      archiveArtifacts(artifacts: 'Testing/Temporary/LastTest.log')
	      throw exc
	    }
	    sh "make install"
	  }
	}
      }
    }
  } }
}

def osxPlatforms = [
  ["gcc", ['CC=gcc-7', 'CXX=g++-7']],
  ["clang", ['CC=/usr/local/opt/llvm/bin/clang', 'CXX=/usr/local/opt/llvm/bin/clang++', 'CXXFLAGS=-I/usr/local/opt/llvm/include', 'LDFLAGS=-L/usr/local/opt/llvm/lib']]
]
for (int i = 0; i < osxPlatforms.size(); i++) {
  def platformEnv = osxPlatforms[i]
  def platform = platformEnv[0]
  platforms["osx-$platform"] = { -> stage("osx-$platform") {
    timeout(time: 1, unit: 'HOURS') {
      node('osx && triqs') {
	def srcDir = pwd()
	def tmpDir = pwd(tmp:true)
	def buildDir = "$tmpDir/build"
	def installDir = "$tmpDir/install"

	dir(installDir) {
	  deleteDir()
	}

	copyArtifacts(projectName: triqsProject, selector: upstream(fallbackToLastSuccessful: true), filter: "osx-${platform}.zip")
	unzip(zipFile: "osx-${platform}.zip", dir: installDir)
	/* fixup zip-stripped permissions (JENKINS-13128) */
	sh "chmod +x $installDir/bin/*"

	checkout scm

	dir(buildDir) { withEnv(platformEnv[1]+[
	    "PATH=$installDir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin",
	    "CPATH=$installDir/include",
	    "LIBRARY_PATH=$installDir/lib",
	    "CMAKE_PREFIX_PATH=$installDir/share/cmake"]) {
	  deleteDir()
	  sh "cmake $srcDir -DTRIQS_ROOT=$installDir"
	  sh "make -j2"
	  try {
	    sh "make test"
	  } catch (exc) {
	    archiveArtifacts(artifacts: 'Testing/Temporary/LastTest.log')
	    throw exc
	  }
	  sh "make install"
	} }
	// zip(zipFile: "osx-${platform}.zip", archive: true, dir: installDir)
      }
    }
  } }
}

parallel platforms
